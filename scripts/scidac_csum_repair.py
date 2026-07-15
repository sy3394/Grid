#!/usr/bin/env python3
"""Check and repair SciDAC checksums in Grid-written LIME files.

Grid builds whose bare crc32 symbol was shadowed by a non-zlib
implementation (oneAPI/SYCL link on Aurora) wrote wrong suma/sumb values
into scidac-checksum records.  The field payload itself is correct.  This
tool recomputes the checksums from the ildg-binary-data / scidac-binary-data
payload (zlib crc32 per site, rotated by global_site%29 / %31, XOR-combined;
see Grid/parallelIO/BinaryIO.h ScidacChecksum) and patches the hex fields
inside the scidac-checksum XML records.  The payload is never modified.

Usage:
    scidac_csum_repair.py file.scidac ...            # check only
    scidac_csum_repair.py --fix file.scidac ...      # patch bad checksums
    scidac_csum_repair.py --fix --backup file ...    # keep file.orig

Lattice geometry and site size are read from the scidac-private-file-xml
and scidac-private-record-xml records, so the tool works for any lattice
volume and any per-site type.  Exit status: 0 all checksums good (or
fixed), 1 mismatches found (check mode) or unfixable, 2 parse error.
"""
import argparse
import os
import re
import struct
import sys
import zlib

LIME_MAGIC = 0x456789AB
LIME_HDR = 144  # 4 magic + 2 version + 2 flags + 8 nbytes + 128 type

BINARY_TYPES = ("ildg-binary-data", "scidac-binary-data")


def lime_records(blob):
    """Yield (hdr_off, rec_type, data_off, nbytes) for each LIME record."""
    off = 0
    while off + LIME_HDR <= len(blob):
        magic, _ver, _flags = struct.unpack_from(">IHH", blob, off)
        if magic != LIME_MAGIC:
            raise ValueError("bad LIME magic %#x at offset %d" % (magic, off))
        (nbytes,) = struct.unpack_from(">Q", blob, off + 8)
        rec_type = blob[off + 16:off + LIME_HDR].split(b"\0")[0].decode()
        yield off, rec_type, off + LIME_HDR, nbytes
        off += LIME_HDR + ((nbytes + 7) // 8) * 8
    if off != len(blob):
        raise ValueError("trailing garbage after offset %d" % off)


def xml_field(xml, tag):
    m = re.search(("<%s>(.*?)</%s>" % (tag, tag)).encode(), xml, re.S)
    return m.group(1).strip().decode() if m else None


def scidac_checksum(payload, sitesize):
    """suma,sumb per QIO/SciDAC: per-site crc32, rotate-left by site%29/%31, XOR."""
    suma = sumb = 0
    nsites = len(payload) // sitesize
    for g in range(nsites):
        crc = zlib.crc32(payload[g * sitesize:(g + 1) * sitesize]) & 0xFFFFFFFF
        r29, r31 = g % 29, g % 31
        suma ^= ((crc << r29) | (crc >> (32 - r29))) & 0xFFFFFFFF
        sumb ^= ((crc << r31) | (crc >> (32 - r31))) & 0xFFFFFFFF
    return suma, sumb


def patch_hex_field(xml, tag, value):
    """Replace <tag>hex</tag> in xml bytes, preserving field width if possible.

    Returns new xml bytes.  Pads with leading zeros up to the old width
    (Grid parses with stoull(...,16), so leading zeros are harmless); grows
    the field only if the new value needs more digits than the old one had.
    """
    pat = ("<%s>([0-9a-fA-F]+)</%s>" % (tag, tag)).encode()
    m = re.search(pat, xml)
    if not m:
        raise ValueError("no <%s> hex field in checksum XML" % tag)
    new = "%x" % value
    width = max(len(m.group(1)), len(new))
    rep = ("<%s>%s</%s>" % (tag, new.zfill(width), tag)).encode()
    return xml[:m.start()] + rep + xml[m.end():]


def rebuild_record(blob, hdr_off, nbytes_old, new_data):
    """Return blob with the record at hdr_off carrying new_data instead."""
    hdr = bytearray(blob[hdr_off:hdr_off + LIME_HDR])
    struct.pack_into(">Q", hdr, 8, len(new_data))
    pad_new = (-len(new_data)) % 8
    old_end = hdr_off + LIME_HDR + ((nbytes_old + 7) // 8) * 8
    return (blob[:hdr_off] + bytes(hdr) + new_data + b"\0" * pad_new
            + blob[old_end:])


def process(path, fix, backup):
    with open(path, "rb") as f:
        blob = f.read()

    dims = None
    sitesize = None
    payload_span = None       # (data_off, nbytes) of last binary-data record
    bad = []                  # (hdr_off, data_off, nbytes, new_xml)
    n_checked = 0

    for hdr_off, rec_type, data_off, nbytes in lime_records(blob):
        data = blob[data_off:data_off + nbytes]
        if rec_type == "scidac-private-file-xml":
            dims = [int(x) for x in xml_field(data, "dims").split()]
        elif rec_type == "grid-format":
            # FieldMetaData: <dimension><elem>N</elem>...</dimension>
            dimension = xml_field(data, "dimension")
            if dimension is not None:
                dims = [int(x) for x in re.findall(r"<elem>(\d+)</elem>", dimension)]
        elif rec_type == "scidac-private-record-xml":
            sitesize = int(xml_field(data, "typesize")) * int(xml_field(data, "datacount"))
        elif rec_type in BINARY_TYPES:
            payload_span = (data_off, nbytes)
        elif rec_type == "scidac-checksum":
            n_checked += 1
            if payload_span is None or sitesize is None or dims is None:
                print("%s: checksum record #%d lacks preceding geometry/payload; skipped"
                      % (path, n_checked))
                continue
            p_off, p_bytes = payload_span
            vol = 1
            for d in dims:
                vol *= d
            if vol * sitesize != p_bytes:
                print("%s: payload size %d != dims %s x sitesize %d; skipped"
                      % (path, p_bytes, dims, sitesize))
                continue
            suma, sumb = scidac_checksum(blob[p_off:p_off + p_bytes], sitesize)
            stored_a = int(xml_field(data, "suma"), 16)
            stored_b = int(xml_field(data, "sumb"), 16)
            if (suma, sumb) == (stored_a, stored_b):
                print("%s: checksum #%d OK (suma=%08x sumb=%08x)"
                      % (path, n_checked, suma, sumb))
            else:
                print("%s: checksum #%d MISMATCH stored suma=%08x sumb=%08x, "
                      "recomputed suma=%08x sumb=%08x"
                      % (path, n_checked, stored_a, stored_b, suma, sumb))
                new_xml = patch_hex_field(data, "suma", suma)
                new_xml = patch_hex_field(new_xml, "sumb", sumb)
                bad.append((hdr_off, data_off, nbytes, new_xml))
            payload_span = None

    if not bad:
        return 0
    if not fix:
        return 1

    if backup:
        with open(path + ".orig", "wb") as f:
            f.write(blob)

    grown = [b for b in bad if len(b[3]) != b[2]]
    if not grown:
        # same-width XML: patch the checksum records in place
        with open(path, "r+b") as f:
            for _hdr_off, data_off, nbytes, new_xml in bad:
                f.seek(data_off)
                f.write(new_xml)
    else:
        # a hex field grew: rebuild the affected records (back to front so
        # earlier offsets stay valid), then rewrite the file atomically
        for hdr_off, _data_off, nbytes, new_xml in sorted(bad, reverse=True):
            blob = rebuild_record(blob, hdr_off, nbytes, new_xml)
        tmp = path + ".tmp"
        with open(tmp, "wb") as f:
            f.write(blob)
        os.replace(tmp, path)
    print("%s: fixed %d checksum record(s)" % (path, len(bad)))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fix", action="store_true", help="patch bad checksums in place")
    ap.add_argument("--backup", action="store_true", help="with --fix, keep <file>.orig")
    ap.add_argument("files", nargs="+")
    args = ap.parse_args()

    status = 0
    for path in args.files:
        try:
            status = max(status, process(path, args.fix, args.backup))
        except (ValueError, OSError) as e:
            print("%s: ERROR %s" % (path, e))
            status = 2
    return status


if __name__ == "__main__":
    sys.exit(main())

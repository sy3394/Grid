"""
Convert qlat topo-field pickle → SCIDAC/LIME files readable by Grid's ScidacReader.

qlat pickle stores fields in (x,y,z,t) order (x slowest, t fastest in memory).
Grid's lexicographic order is x fastest, t slowest.
The conversion is simply .T (reverses all 4 axes before C-order flatten).

Record structure matches Grid's ScidacWriter::writeScidacFieldRecord:
  1. grid-format              (MB=1,ME=0)  FieldMetaData XML
  2. scidac-record-xml        (MB=0,ME=0)  emptyUserRecord XML
  3. scidac-private-record-xml(MB=0,ME=0)  scidacRecord XML
  4. ildg-binary-data         (MB=0,ME=0)  big-endian complex128 payload
  5. scidac-checksum          (MB=0,ME=1)  SciDAC checksum XML
"""

import struct
import zlib
import pickle
import numpy as np
import sys
import os

# ──────────────────────────────────────────────────────────────────────────────
# LIME record writer
# ──────────────────────────────────────────────────────────────────────────────

LIME_MAGIC   = 0x456789AB   # 32-bit magic (first 4 bytes of the header word)
LIME_VERSION = 1

def _lime_record(type_str: str, data: bytes, MB: bool, ME: bool) -> bytes:
    """Pack one LIME record (144-byte header + data padded to 8-byte boundary).

    Binary header layout (from lime_binary_header.h, 18 × 8-byte words = 144 bytes):
      bytes  0- 3 : magic       = 0x456789AB  (uint32)
      bytes  4- 5 : version     = 1           (uint16)
      byte   6    : MB|ME flags               (uint8,  0x80=MB 0x40=ME)
      byte   7    : reserved    = 0           (uint8)
      bytes  8-15 : data_length              (uint64)
      bytes 16-143: record type string        (128 bytes, null-padded)
    """
    flags_byte = (0x80 if MB else 0) | (0x40 if ME else 0)
    type_b     = type_str.encode('ascii').ljust(128, b'\x00')[:128]
    header     = struct.pack('>IHBBQ128s', LIME_MAGIC, LIME_VERSION,
                             flags_byte, 0, len(data), type_b)
    # 4+2+1+1+8+128 = 144 bytes
    pad_len    = (8 - len(data) % 8) % 8
    return header + data + b'\x00' * pad_len


# ──────────────────────────────────────────────────────────────────────────────
# SciDAC checksum (mirrors BinaryIO.h::ScidacChecksum)
# ──────────────────────────────────────────────────────────────────────────────

def _rotl32(x: int, n: int) -> int:
    """Rotate 32-bit integer x left by n bits."""
    n &= 31
    return ((x << n) | (x >> (32 - n))) & 0xFFFFFFFF if n else (x & 0xFFFFFFFF)


def _scidac_checksums(payload: bytes, Nsites: int):
    """
    Compute Grid's SciDAC checksums over big-endian binary field data.

    Grid calls ScidacChecksum *after* htobe (on WRITE) and *before* be64toh
    (on READ), so the CRC is taken over the big-endian bytes on disk — which
    is exactly what we have in `payload`.

    Sites must be in Grid's lexicographic order: x fastest, t slowest.
    Returns (csuma, csumb) as unsigned 32-bit integers.
    """
    mv   = memoryview(payload)    # zero-copy slicing
    csuma = csumb = 0
    for site in range(Nsites):
        crc   = zlib.crc32(mv[site * 16 : site * 16 + 16]) & 0xFFFFFFFF
        csuma ^= _rotl32(crc, site % 29)
        csumb ^= _rotl32(crc, site % 31)
    return csuma, csumb


# ──────────────────────────────────────────────────────────────────────────────
# Grid-compatible SCIDAC writer for a complex scalar field
# ──────────────────────────────────────────────────────────────────────────────

def write_scidac_complex(fname: str, field: np.ndarray):
    """
    Write a 4D complex scalar field to a Grid-compatible SCIDAC/LIME file.

    Parameters
    ----------
    fname : output filename
    field : numpy array, shape (Nt, Nz, Ny, Nx), dtype complex128
            Must be in Grid order: t slowest, x fastest.
            (Pass arr.T when arr comes from qlat's (x,y,z,t) layout.)
    """
    assert field.ndim == 4, "field must be 4D"
    Nt, Nz, Ny, Nx = field.shape
    Nsites = Nx * Ny * Nz * Nt

    # ── Binary payload: big-endian interleaved re/im ──────────────────────────
    # Ravel in C-order over (Nt,Nz,Ny,Nx) → t slowest, x fastest = Grid order
    flat    = np.ascontiguousarray(field, dtype=np.complex128).ravel()
    be      = np.empty(2 * len(flat), dtype='>f8')
    be[0::2] = flat.real
    be[1::2] = flat.imag
    payload = be.tobytes()          # 16 * Nsites bytes

    # ── SciDAC checksums ──────────────────────────────────────────────────────
    print(f'  computing SciDAC checksums over {Nsites:,} sites…', end=' ', flush=True)
    csuma, csumb = _scidac_checksums(payload, Nsites)
    print(f'done  csuma=0x{csuma:08x}  csumb=0x{csumb:08x}')

    # ── XML records ──────────────────────────────────────────────────────────
    # Record 1: grid-format — FieldMetaData (Grid's proprietary header)
    # Field order in XML must match GRID_SERIALIZABLE_CLASS_MEMBERS in MetaData.h
    field_meta_xml = (
        '<?xml version="1.0"?>'
        '<FieldMetaData>'
          '<nd>4</nd>'
          f'<dimension>'
            f'<elem>{Nx}</elem><elem>{Ny}</elem><elem>{Nz}</elem><elem>{Nt}</elem>'
          f'</dimension>'
          '<boundary>'
            '<elem>PERIODIC</elem><elem>PERIODIC</elem>'
            '<elem>PERIODIC</elem><elem>PERIODIC</elem>'
          '</boundary>'
          '<data_start>0</data_start>'
          '<hdr_version></hdr_version>'
          '<storage_format></storage_format>'
          '<link_trace>0</link_trace>'
          '<plaquette>0</plaquette>'
          '<checksum>0</checksum>'
          f'<scidac_checksuma>{csuma}</scidac_checksuma>'
          f'<scidac_checksumb>{csumb}</scidac_checksumb>'
          '<sequence_number>0</sequence_number>'
          '<data_type>GRID_D_Complex</data_type>'
          '<ensemble_id></ensemble_id>'
          '<ensemble_label></ensemble_label>'
          '<ildg_lfn></ildg_lfn>'
          '<creator>pickle_to_scidac.py</creator>'
          '<creator_hardware></creator_hardware>'
          '<creation_date></creation_date>'
          '<archive_date></archive_date>'
          '<floating_point>IEEE64BIG</floating_point>'
        '</FieldMetaData>'
    ).encode()

    # Record 2: scidac-record-xml — emptyUserRecord
    user_xml = (
        '<?xml version="1.0"?>'
        '<emptyUserRecord><dummy>0</dummy></emptyUserRecord>'
    ).encode()

    # Record 3: scidac-private-record-xml — scidacRecord
    # datatype "GRID_D_Complex" matches ScidacRecordTypeString<iScalar<iScalar<iScalar<ComplexD>>>>
    scidac_rec_xml = (
        '<?xml version="1.0"?>'
        '<scidacRecord>'
          '<version>1</version>'
          '<date></date>'
          '<recordtype>0</recordtype>'
          '<datatype>GRID_D_Complex</datatype>'
          '<precision>D</precision>'
          '<colors>1</colors>'
          '<spins>1</spins>'
          '<typesize>16</typesize>'
          '<datacount>1</datacount>'
        '</scidacRecord>'
    ).encode()

    # Record 5: scidac-checksum — scidacChecksum (hex strings, no "0x" prefix)
    checksum_xml = (
        '<?xml version="1.0"?>'
        '<scidacChecksum>'
          '<version>1</version>'
          f'<suma>{csuma:x}</suma>'
          f'<sumb>{csumb:x}</sumb>'
        '</scidacChecksum>'
    ).encode()

    # ── Write LIME file ───────────────────────────────────────────────────────
    with open(fname, 'wb') as f:
        f.write(_lime_record('grid-format',                 field_meta_xml, MB=True,  ME=False))
        f.write(_lime_record('scidac-record-xml',           user_xml,       MB=False, ME=False))
        f.write(_lime_record('scidac-private-record-xml',   scidac_rec_xml, MB=False, ME=False))
        f.write(_lime_record('ildg-binary-data',            payload,        MB=False, ME=False))
        f.write(_lime_record('scidac-checksum',             checksum_xml,   MB=False, ME=True))

    size_mb = os.path.getsize(fname) / 1024 / 1024
    print(f'  wrote {fname}  ({size_mb:.1f} MB)')


# ──────────────────────────────────────────────────────────────────────────────
# Main: convert all topo_field_* entries in the pickle
# ──────────────────────────────────────────────────────────────────────────────

def main():
    pickle_path = sys.argv[1] if len(sys.argv) > 1 else '32c-hmc-test-demo.pickle'
    out_dir     = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(pickle_path)

    print(f'Loading {pickle_path} …')
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)

    for key, arr in data.items():
        print(f'\n{key}: shape={arr.shape}  dtype={arr.dtype}')
        print(f'  Q = {arr.real.sum():.6f}  (imag rms = {np.sqrt(np.mean(arr.imag**2)):.3e})')

        # qlat: (x,y,z,t) x-slowest, t-fastest in memory (C-order)
        # Grid: x-fastest, t-slowest → .T reverses all axes
        field_grid = arr.T                  # (Nx,Ny,Nz,Nt) → (Nt,Nz,Ny,Nx)
        assert field_grid.ndim == 4

        out_name = os.path.join(out_dir, f'{key}.scidac')
        write_scidac_complex(out_name, field_grid)


if __name__ == '__main__':
    main()

"""
Convert qlat topo-field pickle → SCIDAC/LIME files readable by Grid's ScidacReader.

qlat pickle stores fields in (x,y,z,t) order (x slowest).
Grid/SCIDAC expects (t,z,y,x) order (t slowest, x fastest) in the binary payload.
The conversion is simply .T (reverses all 4 axes).
"""

import struct
import pickle
import numpy as np
import sys
import os

# ──────────────────────────────────────────────────────────────────────────────
# LIME record writer
# ──────────────────────────────────────────────────────────────────────────────

LIME_MAGIC   = 0x456789AB01234567
LIME_VERSION = 1

def _lime_record(type_str: str, data: bytes, MB: bool, ME: bool) -> bytes:
    """Pack one LIME record (header + data padded to 8-byte boundary)."""
    flags     = (0x8000 if MB else 0) | (0x4000 if ME else 0)
    type_b    = type_str.encode('ascii').ljust(128, b'\x00')[:128]
    header    = struct.pack('>QHHQ128s', LIME_MAGIC, LIME_VERSION,
                            flags, len(data), type_b)
    pad_len   = (8 - len(data) % 8) % 8
    return header + data + b'\x00' * pad_len


# ──────────────────────────────────────────────────────────────────────────────
# SCIDAC writer for a single complex scalar field
# ──────────────────────────────────────────────────────────────────────────────

def write_scidac_complex(fname: str, field: np.ndarray):
    """
    Write a 4D complex scalar field to a SCIDAC/LIME file.

    Parameters
    ----------
    fname : output filename
    field : numpy array, shape (Nt, Nz, Ny, Nx), dtype complex128
            Must already be in Grid order: t slowest, x fastest.
    """
    assert field.ndim == 4, "field must be 4D"
    Nt, Nz, Ny, Nx = field.shape

    # ── Record 1: scidac-file-xml  (standalone message) ──────────────────────
    file_xml = (
        '<?xml version="1.0"?>'
        '<scidacFile>'
          '<version>1.1</version>'
          '<spacetime>4</spacetime>'
          f'<dims>{Nx} {Ny} {Nz} {Nt}</dims>'
          '<volfmt>0</volfmt>'
        '</scidacFile>'
    ).encode()
    rec1 = _lime_record('scidac-file-xml', file_xml, MB=True, ME=True)

    # ── Record 2: scidac-record-xml  (begins the field message) ──────────────
    rec_xml = (
        '<?xml version="1.0"?>'
        '<scidacRecord>'
          '<version>1.0</version>'
          '<date></date>'
          '<globaldata>0</globaldata>'
          '<datatype>4D_COMPLEX_DOUBLE</datatype>'
          '<precision>D</precision>'
          '<colors>1</colors>'
          '<spins>1</spins>'
          '<typesize>16</typesize>'   # 2 × float64
          '<datacount>1</datacount>'
        '</scidacRecord>'
    ).encode()
    rec2 = _lime_record('scidac-record-xml', rec_xml, MB=True, ME=False)

    # ── Record 3: scidac-binary-data  (ends the field message) ───────────────
    # Interleave re/im and byte-swap to big-endian
    flat = np.ascontiguousarray(field, dtype=np.complex128).ravel()
    be   = np.empty(2 * len(flat), dtype='>f8')
    be[0::2] = flat.real
    be[1::2] = flat.imag
    rec3 = _lime_record('scidac-binary-data', be.tobytes(), MB=False, ME=True)

    with open(fname, 'wb') as f:
        f.write(rec1)
        f.write(rec2)
        f.write(rec3)

    size_kb = os.path.getsize(fname) / 1024
    print(f'  wrote {fname}  ({size_kb:.0f} kB)')


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

        # qlat stores (x,y,z,t); Grid wants (t,z,y,x) — .T reverses all axes
        field_grid = arr.T                      # (Nx,Ny,Nz,Nt) → (Nt,Nz,Ny,Nx)
        assert field_grid.ndim == 4

        out_name = os.path.join(out_dir, f'{key}.scidac')
        write_scidac_complex(out_name, field_grid)


if __name__ == '__main__':
    main()

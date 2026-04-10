"""
Convert qlat topo-field pickle → raw big-endian binary files readable by
Grid's BinaryIO::readLatticeObject (no LIME wrapper, no HDF5 dependency).

Axis convention
---------------
qlat stores fields as numpy arrays with shape (Nx, Ny, Nz, Nt), x slowest
in C-order (i.e. t varies fastest in memory).
Grid's lexicographic ordering has x fastest, t slowest.
After .T the shape becomes (Nt, Nz, Ny, Nx); C-order ravel then gives
Grid lex order: x varies fastest.

File format
-----------
Raw big-endian IEEE754 double-precision complex values, no header.
Each site is 16 bytes: 8 bytes real (big-endian float64) + 8 bytes imag.
Total size = Nsites * 16 bytes.

On the C++ side Grid reads this with:
    BinarySimpleMunger<ComplexD, ComplexD> munge;
    BinaryIO::readLatticeObject<vobj, ComplexD>(out, fname, munge,
        /*offset=*/0, "IEEE64BIG", nersc_csum, csuma, csumb);

Usage
-----
    python3 pickle_to_binary.py <pickle_file> <output_dir>

Example
-------
    python3 pickle_to_binary.py ../../data/32c-hmc-test-demo.pickle data/
    → data/topo_field_0.bin  data/topo_field_1.bin  (16 MB each)
"""

import pickle
import numpy as np
import sys
import os


def write_binary_complex(fname: str, field: np.ndarray):
    """
    Write a 4D complex scalar field as raw big-endian binary.

    Parameters
    ----------
    fname  : output filename (.bin)
    field  : numpy array, shape (Nt, Nz, Ny, Nx), dtype complex128
             Must be in Grid order: t slowest, x fastest.
             (Pass arr.T when arr comes from qlat's (x,y,z,t) layout.)
    """
    assert field.ndim == 4, "field must be 4D"
    Nt, Nz, Ny, Nx = field.shape
    Nsites = Nx * Ny * Nz * Nt

    # C-order ravel: t slowest, x fastest = Grid lexicographic order.
    flat = np.ascontiguousarray(field, dtype=np.complex128).ravel()

    # Interleave real and imag as big-endian float64.
    # Grid's BinaryIO expects IEEE64BIG: re0, im0, re1, im1, ...
    out = np.empty(2 * Nsites, dtype='>f8')
    out[0::2] = flat.real
    out[1::2] = flat.imag

    with open(fname, 'wb') as f:
        out.tofile(f)

    size_mb = os.path.getsize(fname) / 1024 / 1024
    print(f'  wrote {fname}  ({size_mb:.1f} MB,  {Nsites} sites)')


def main():
    pickle_path = sys.argv[1] if len(sys.argv) > 1 else '32c-hmc-test-demo.pickle'
    out_dir     = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(pickle_path)

    print(f'Loading {pickle_path} ...')
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)

    for key, arr in data.items():
        print(f'\n{key}: shape={arr.shape}  dtype={arr.dtype}')
        print(f'  Q = {arr.real.sum():.6f}  (imag rms = {np.sqrt(np.mean(arr.imag**2)):.3e})')

        # qlat (x,y,z,t) x-slowest → Grid (t,z,y,x) t-slowest via .T
        field_grid = arr.T
        assert field_grid.ndim == 4

        out_name = os.path.join(out_dir, f'{key}.bin')
        write_binary_complex(out_name, field_grid)


if __name__ == '__main__':
    main()

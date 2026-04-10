"""
Convert qlat topo-field pickle → HDF5 files readable by Grid via H5Cpp.

Axis convention
---------------
qlat stores fields as numpy arrays with shape (Nx, Ny, Nz, Nt), x slowest.
Grid's lexicographic ordering has x fastest, t slowest.
After .T (shape → (Nt, Nz, Ny, Nx)), C-order ravel gives Grid order.

HDF5 layout
-----------
Dataset  "field"    shape (Nsites, 2)  dtype float64
                    axis-0: site in Grid lexicographic order (x fastest)
                    axis-1: [0]=real, [1]=imag
Attributes
    Nx, Ny, Nz, Nt  — lattice dimensions (for validation in C++)
"""

import h5py
import pickle
import numpy as np
import sys
import os


def write_hdf5_complex(fname: str, field: np.ndarray):
    """
    Write a 4D complex scalar field to HDF5 in Grid's lexicographic order.

    Parameters
    ----------
    fname  : output filename (.h5)
    field  : numpy array, shape (Nt, Nz, Ny, Nx), dtype complex128
             Grid order: t slowest, x fastest.
             (Pass arr.T when arr comes from qlat's (x,y,z,t) layout.)
    """
    assert field.ndim == 4
    Nt, Nz, Ny, Nx = field.shape
    Nsites = Nx * Ny * Nz * Nt

    # C-order ravel: t slowest, x fastest = Grid lexicographic order
    flat = np.ascontiguousarray(field, dtype=np.complex128).ravel()

    # Interleaved real/imag as float64 rows — shape (Nsites, 2)
    data = np.empty((Nsites, 2), dtype=np.float64)
    data[:, 0] = flat.real
    data[:, 1] = flat.imag

    with h5py.File(fname, 'w') as f:
        ds = f.create_dataset('field', data=data, dtype='float64')
        f.attrs['Nx'] = Nx
        f.attrs['Ny'] = Ny
        f.attrs['Nz'] = Nz
        f.attrs['Nt'] = Nt

    size_mb = os.path.getsize(fname) / 1024 / 1024
    print(f'  wrote {fname}  ({size_mb:.1f} MB)')


def main():
    pickle_path = sys.argv[1] if len(sys.argv) > 1 else '32c-hmc-test-demo.pickle'
    out_dir     = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(pickle_path)

    print(f'Loading {pickle_path} …')
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)

    for key, arr in data.items():
        print(f'\n{key}: shape={arr.shape}  dtype={arr.dtype}')
        print(f'  Q = {arr.real.sum():.6f}  (imag rms = {np.sqrt(np.mean(arr.imag**2)):.3e})')

        # qlat (x,y,z,t) x-slowest → Grid (t,z,y,x) t-slowest via .T
        field_grid = arr.T
        assert field_grid.ndim == 4

        out_name = os.path.join(out_dir, f'{key}.h5')
        write_hdf5_complex(out_name, field_grid)


if __name__ == '__main__':
    main()

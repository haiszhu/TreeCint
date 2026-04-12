#!/usr/bin/env python3
import sys
import numpy as np
import h5py
from scipy.io import savemat
from pyscf import gto


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: verify_eri_h2o_ref.py <eri_h5_filename> <out_mat_filename>")
        return 2

    eri_h5_filename = sys.argv[1]
    out_mat_filename = sys.argv[2]

    geom = '''
    O  0    0.       0.
    H  0    -0.757   0.587
    H  0    0.757    0.587
    '''

    mol = gto.M(verbose=0, atom=geom, basis='ccpvdz')
    ref = mol.intor('int2e')

    with h5py.File(eri_h5_filename, 'r') as f:
        eri = f['DS1'][()]

    d = eri - ref
    max_abs = float(np.max(np.abs(d)))
    rel_fro = float(np.linalg.norm(d.ravel()) / np.linalg.norm(ref.ravel()))
    rmse = float(np.sqrt(np.mean(d * d)))

    savemat(out_mat_filename, {
        'max_abs': np.array([[max_abs]], dtype=np.float64),
        'rel_fro': np.array([[rel_fro]], dtype=np.float64),
        'rmse': np.array([[rmse]], dtype=np.float64),
    })
    print(f"saved metrics to {out_mat_filename}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

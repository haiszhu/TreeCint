# TreeCint

A MATLAB-based pipeline for computing **Electron Repulsion Integrals (ERIs)** via the
**Interpolative Separable Density Fitting (ISDF)** method, accelerated by an adaptive
octree Coulomb solver (Box-DMK).

## Provenance

The original version of this code path lives in
[`haiszhu/libcint`](https://github.com/haiszhu/libcint/tree/master), a fork and
extension of the upstream
[`sunqm/libcint`](https://github.com/sunqm/libcint) library.

This repository modernizes that original workflow into a more reproducible
MATLAB-based pipeline for ERI construction via ISDF and Box-DMK with support from Codex and Claude Code.

The modernization introduces reproducible external dependencies and clearer build/integration
boundaries, including:
- `external/libcint` (GTO evaluation backend) - https://github.com/sunqm/libcint
- `external/libid` (interpolative decomposition backend, legacy/default path) - https://github.com/fastalgorithms/libid
- `external/coqui` (interpolative decomposition backend, planned replacement for `libid`) - https://github.com/AbInitioQHub/coqui
- `external/dmk` (Box-DMK Coulomb kernel backend) - https://github.com/flatironinstitute/dmk
- `external/treefun` (adaptive MATLAB octree representation) - https://github.com/haiszhu/treefun

## Paper Context

This implementation supports: Zhu, H., Yeh, C.-N., Morales, M. A., Greengard, L., Jiang, S., and Kaye, J. *Interpolative separable density fitting on adaptive real space grids*. arXiv:2510.20826 (2026). [https://arxiv.org/pdf/2510.20826](https://arxiv.org/pdf/2510.20826)

The pipeline avoids explicit O(N⁴) four-index ERI storage by:
1. Resolving orbital products on an adaptive octree (treefun3)
2. Compressing via interpolative decomposition (ISDF / libid)
3. Evaluating the Coulomb kernel on the compressed representation (Box-DMK)
4. Contracting back to the full four-index tensor

## Repository Layout

```
TreeCint-legacy/
├── Makefile                        # Top-level build (all targets)
├── basis/
│   └── cc-pvdz.dat                 # NWChem-format cc-pVDZ basis set
├── external/
│   ├── libcint/                    # GTO integral library (C, CMake submodule)
│   ├── libid/                      # Interpolative decomposition (Fortran submodule)
│   ├── dmk/                        # Box-DMK Coulomb solver (Fortran submodule)
│   └── treefun/                    # Adaptive octree in MATLAB (submodule)
├── src/
│   ├── wrapper/
│   │   ├── gateway.mw              # MWrap spec → TreeCint.mex (GTO evaluator)
│   │   ├── cgto.c                  # C bridge between MWrap and libcint
│   │   └── cgto.h
│   ├── bdmk/
│   │   ├── bdmk.mw                 # MWrap spec → bdmk_mex.mex (Coulomb solver)
│   │   └── bdmk_mex.f90            # Fortran wrapper looping over nd density components
│   └── libid/
│       ├── libid.mw                # MWrap spec → libid_mex.mex (ID solver)
│       └── libid_mex.f90           # Fortran bridge to iddr_aid
├── matlab/                         # MATLAB utilities + compiled .mexmaca64 binaries
│   ├── gto.m                       # MATLAB molecule builder (mirrors PySCF gto.M)
│   ├── treefun2bdmk.m              # Convert treefun3 octree → BDMK format
│   ├── computeVmunu_bdmk.m         # Assemble Coulomb kernel matrix via bdmk_mex
│   ├── computeVijkl.m              # Contract Vmunu → four-index ERI tensor
│   ├── id_libid.m                  # Rank-k interpolative decomposition wrapper
│   ├── loadNWChemBasis.m           # NWChem basis file parser
│   ├── write_treefun_h5.m          # Write tree metadata to HDF5
│   ├── TreeCint.mexmaca64          # Compiled GTO evaluator (Apple Silicon)
│   ├── bdmk_mex.mexmaca64          # Compiled BDMK Coulomb solver (Apple Silicon)
│   └── libid_mex.mexmaca64         # Compiled ID solver (Apple Silicon)
└── test/
    └── isdf/
        └── H2O/
            └── test_ccpvdz.m       # Canonical H2O/cc-pVDZ end-to-end test
```

## External Dependencies

| Submodule | Language | Build | Role |
|-----------|----------|-------|------|
| `external/libcint` | C | CMake → `libcint.dylib` | Evaluate spherical GTO basis functions at grid points |
| `external/libid` | Fortran | Make → `id_lib.a` | Randomized interpolative decomposition (`iddr_aid`) |
| `external/dmk` | Fortran | Make → `libbdmk_ref.a` | Box-DMK: tree-based Coulomb kernel `1/|r−r'|` |
| `external/treefun` | MATLAB | — | Adaptive octree function representation (`treefun3`) |

## Build Requirements

- **macOS (Apple Silicon / maca64)**
- `gcc-14` (C compiler)
- `gfortran` (Fortran compiler)
- `cmake`
- `mwrap` at `~/mwrap/mwrap` (override with `MWRAP=...`)
- MATLAB (auto-detected from `/Applications/MATLAB_R*.app`)
- Python with `numpy`, `h5py`, `pyscf`, `scipy` (for verification step)

## Build Targets

```bash
make libcint-all       # CMake configure + build + install libcint
make libid-build       # Build libid static archive and run upstream tests
make libid-mex         # Build libid_mex.mexmaca64 via MWrap
make bdmk-lib          # Compile 22 Fortran sources → libbdmk_ref.a
make bdmk-mex          # Build bdmk_mex.mexmaca64 via MWrap
make treecint-mex      # Build TreeCint.mexmaca64 (depends on libcint + libid)
make all               # treecint-mex + bdmk-mex (full default build)

make clean             # Remove all build artifacts
```

### OpenMP

```bash
make all USE_OPENMP=1   # Enable OpenMP in GTO evaluation and BDMK
```

## How MWrap Works

Each `.mw` file is the source of truth for a MEX interface. MWrap is run twice per target:

```bash
mwrap -c99complex -mex TreeCint -mb -list gateway.mw   # → MATLAB stub .m files
mwrap -c99complex -mex TreeCint -c TreeCint_gateway.c gateway.mw  # → MEX C glue
```

The MATLAB stubs (e.g. `GTOval_sph_generic_mwrap_mex.m`) handle array-size bookkeeping
before calling into the compiled MEX binary.

## Pipeline Overview

The full pipeline is demonstrated in `test/isdf/H2O/test_ccpvdz.m` for H₂O / cc-pVDZ
(Norb = 24).

```
H₂O geometry + cc-pVDZ basis
        │
        ▼  gto.m  (MATLAB GTO builder, mirrors PySCF gto.M)
   mol struct  (atm, bas, env arrays in PySCF layout)
        │
        ▼  treefun3 + mol.eval_gto → TreeCint.mex (libcint GTOval_sph)
   adaptive octree f  resolving all φ_i(r) to tolerance ε
        │
        ▼  treefun2bdmk  (upsample + Gauss-Legendre quadrature + vol_tree_fix_lr)
   BDMK tree data: src, itree, iptr, centers, boxsize, wtsleaf
        │
        ▼  mol.eval_gto at quadrature points
   fvals_ij[npts × Norb*(Norb+1)/2]   orbital products φ_i(r)·φ_j(r)
        │
        ▼  id_libid → libid_mex.mex (iddr_aid, rank-nd ID)
   Ask[npts × nd]  interpolating vectors
   idcoefs[nd × Norb*(Norb+1)/2]  ISDF coefficients
        │
        ▼  computeVmunu_bdmk → bdmk_eval_mex → bdmk_mex.mex (dmk bdmk.f)
   Vmunu[nd × nd]   2-index Coulomb kernel in ISDF basis
        │
        ▼  computeVijkl  (matrix contraction)
   Vijkl[Norb × Norb × Norb × Norb]   electron repulsion integrals
        │
        ▼  verify_eri_with_python → PySCF int2e reference
   max_abs / rel_fro / rmse
```

### Key Parameters (H₂O / cc-pVDZ test)

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `treefun_order` | 8 | Chebyshev/Legendre order per box dimension (8³ = 512 pts/box) |
| `treefun_eps` | 1e-3 / 1e-6 | Adaptive tree tolerance |
| `isdf_eps` | 1e-3 / 1e-6 | ISDF rank target |
| `nd` | 290 | Number of ISDF interpolation points |
| `rad` | 15 | Box half-width in Bohr |
| `ikernel` | 1 | Coulomb kernel `1/|r−r'|` |
| `beta` | 6.0 | BDMK plane-wave parameter |
| `ipoly` | 0 | Legendre quadrature nodes |

### Accuracy (from test comments)

| Tolerance | max\_abs | rel\_fro | rmse |
|-----------|----------|----------|------|
| 1e-3 | 8.3e-5 | 6.3e-5 | 2.9e-6 |
| 1e-6 | 1.5e-8 | 1.5e-8 | 6.7e-10 |

## Module Details

### `gto.m`

MATLAB reimplementation of PySCF's `gto.M`. Parses XYZ geometry, loads NWChem basis files,
constructs `atm[natm×6]` / `bas[nbas×8]` / `env[...]` arrays in C-compatible (row-major)
layout with full GTO normalization. Returns a `mol` struct exposing:

- `mol.eval_gto(name, coords)` — all Norb orbitals at grid points
- `mol.eval_gto2(name, coords)` — all Norb*(Norb+1)/2 orbital products
- `mol.eval_gtoi(name, coords, i)` — the i-th orbital only

### `cgto.c`

C bridge between MWrap-generated MEX glue and libcint. The generic entry points
(`GTOval_sph_generic_mwrap`, `GTOival_sph_generic_mwrap`) dynamically size the
non-zero screening table (`non0tab`) to `nbas` entries of 1 (no screening) and
evaluate libcint's `GTOval_sph` one grid point at a time inside an OpenMP loop.

### `bdmk_mex.f90`

Fortran wrapper that calls the upstream `bdmk` subroutine in a loop over the `nd`
density components. Each call computes the full tree traversal for one density, then
accumulates the result into `pot(nd, npbox, nboxes)`.

### `computeVmunu_bdmk.m`

Assembles the nd×nd Coulomb kernel matrix:

```
Vmunu[μ,ν] = ∫∫ Ask_μ(r) · (1/|r−r'|) · Ask_ν(r') dr dr'
```

BDMK returns potentials in scaled coordinates; the wrapper applies `/ ratio²` to recover
physical units, then integrates with leaf quadrature weights.

### `computeVijkl.m`

Contracts ISDF representation back to four-index ERIs:

```
Vijkl[i,j,k,l] = Σ_{μ,ν} C_{μ,ij} · Vmunu[μ,ν] · C_{ν,kl}
```

Implemented as a double loop over `(i, k)` with matrix multiplies — O(Norb² · nd²).

## Python Verification

The canonical H2O test runs a Python/PySCF `int2e` reference check and reports max absolute
error, relative Frobenius norm, and RMSE against the generated ERI tensor.

Required packages:
```bash
python -m pip install numpy h5py pyscf scipy
```

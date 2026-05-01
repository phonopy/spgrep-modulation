# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`spgrep-modulation` analyzes collective atomic modulations using irreducible space-group representations of phonon eigenmodes, computes isotropy subgroups on the fly, and generates modulated structures along order-parameter directions. Built on top of `phonopy` (dynamical matrices), `spgrep` (space-group irreps), and `hsnf` (Hermite/Smith normal forms).

## Commands

The project uses `uv` and a `justfile`:

- `just install` — `uv sync --all-extras` (sets up `.venv` with dev/docs/vis extras)
- `just test` — `uv run pytest -v`
- `just prek` — `prek run --all-files` (pre-commit equivalent; run before committing)
- `just docs` — `uv run sphinx-autobuild docs docs_build`

Other:

- Single test: `uv run pytest tests/test_modulation.py::test_name -v`
- Coverage (matches CI): `uv run pytest --cov=src/spgrep_modulation --cov-report=xml tests/`
- Build docs once: `uv run sphinx-build docs docs_build`

CI (`.github/workflows/testing.yml`) runs Python 3.10–3.13 and **installs `spgrep` from its `develop` branch** rather than PyPI. If a test fails locally but passes in CI (or vice versa), suspect a `spgrep` API skew and check `develop` of `spglib/spgrep`.

## Architecture

The package has a single public entry point — `Modulation` — and four supporting modules. The data flow forms a pipeline; understanding the pipeline is the key to navigating the code.

```
phonopy DynamicalMatrix + qpoint
    │
    ▼
irreps.get_eigenmode_representation        # build rep on phonon eigenmodes
    │
    ▼
irreps.project_eigenmode_representation    # decompose into irreps via spgrep
    │
    ▼
modulation.Modulation._group_eigenspaces   # block-diagonalize modified DM per irrep
    │   produces self._eigenspaces: list of (eigval, basis (dim_irrep, num_atoms, 3), irrep_matrices)
    │
    ▼
isotropy.IsotropyEnumerator                # enumerate isotropy subgroups of each irrep
    │
    ▼
modulation.Modulation.get_high_symmetry_modulated_supercells
                                           # combine basis + isotropy directions into PhonopyAtoms
```

Module roles:

- `modulation.py` — `Modulation` class. Construct via `Modulation.with_supercell_and_symmetry_search(dynamical_matrix, supercell_matrix, qpoint, factor=...)`. After init, `eigenspaces[i] = (eigval, basis, irrep_repmat)`. `get_high_symmetry_modulated_supercells(i)` returns `PhonopyAtoms` along 1D order-parameter directions. The phonopy convention swap (`phase = exp(2πi q·r)`) appears around `_group_eigenspaces`; the "modified dynamical matrix" (see `utils.get_modified_dynamical_matrix`) is what makes per-irrep block-diagonalization correct.
- `irreps.py` — Bridges phonopy's eigenmode representation to spgrep's `project_to_irrep`. Returns `(all_basis, irreps, mapping_little_group)`.
- `isotropy.py` — `IsotropyEnumerator` enumerates isotropy subgroups of a small representation at a qpoint. Uses HNF/SNF (`hsnf`) to handle integer linear systems over translations and `spgrep.symmetry.subgroup.enumerate_point_subgroup` for point-group subgroups.
- `utils.py` — Numerical helpers: `qr_unique` (sign-stable QR), `get_modified_dynamical_matrix`, `sample_on_unit_sphere`, `gcd_on_list`/`lcm_on_list`, and the `NDArray*` aliases used as type hints throughout.
- `visualize.py` — Optional (extras `[vis]`); requires `nglview`/`pymatgen`. Don't import from `__init__` paths if working in a minimal env.

### Test fixtures

`tests/conftest.py` exposes session-scoped phonopy fixtures: `ph_bto` (BaTiO3, mp-2998), `ph_mgo` (mp-1265, Fm-3m), `ph_si_diamond` (mp-149, Fd-3m), `ph_aln` (mp-661, P6_3mc). The `.yaml.xz` files come from phonondb (CC BY 4.0) — do not regenerate them.

## Conventions specific to this repo

- Python ≥ 3.10. `from __future__ import annotations` is used in modules with complex type hints.
- Ruff lint: `E`, `F`, `I`, `UP` enabled, `E501` ignored (line length 99). Most `D*` pydocstyle rules are disabled — don't add docstrings beyond what's already there just to satisfy a linter.
- Numpy array shapes are documented inline in comments at non-trivial einsum/reshape sites; preserve and update them when editing those blocks.
- The codebase uses `np.einsum(..., optimize="greedy")` deliberately for the irrep projection hot paths; don't replace with naive matmul without benchmarking.
- `setuptools-scm` provides the version — never hardcode `__version__`.
- `[tool.uv] exclude-newer = "1 week"` pins the resolver to packages released ≥ 1 week ago. If a fresh `uv sync` produces unexpectedly old versions, this is why.

## Working with this repo

- Active branch is typically `develop`; `main` is the release branch. Both are CI-tested.
- The `develop` branch publishes docs to `phonopy.github.io/spgrep-modulation/develop/`.
- Untracked scratch files (`debug.py`, `rep.npy`, ad hoc notebooks under `docs/`) may exist in the working tree — leave them alone unless asked.

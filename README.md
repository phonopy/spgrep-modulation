# spgrep-modulation
[![testing](https://github.com/phonopy/spgrep-modulation/actions/workflows/testing.yml/badge.svg)](https://github.com/phonopy/spgrep-modulation/actions/workflows/testing.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/phonopy/spgrep-modulation/main.svg)](https://results.pre-commit.ci/latest/github/phonopy/spgrep-modulation/main)
[![codecov](https://codecov.io/gh/lan496/spgrep-modulation/branch/main/graph/badge.svg?token=K80FQQJ383)](https://codecov.io/gh/lan496/spgrep-modulation)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spgrep-modulation)](https://img.shields.io/pypi/pyversions/spgrep-modulation)
[![PyPI version](https://badge.fury.io/py/spgrep-modulation.svg)](https://badge.fury.io/py/spgrep-modulation)
[![PyPI Downloads](https://img.shields.io/pypi/dm/spgrep-modulation)](https://img.shields.io/pypi/dm/spgrep-modulation)

Collective atomic modulation analysis with irreducible space-group representation

- Github: <https://github.com/phonopy/spgrep-modulation>
- Document: <https://phonopy.github.io/spgrep-modulation>
- Document (develop): <https://phonopy.github.io/spgrep-modulation/develop/>
- PyPI: <https://pypi.org/project/spgrep-modulation>

## Features

- Calculate representation matrices and irreps formed by phonon eigenmodes
- Calculate isotropy subgroups of irreps of space groups on the fly
- Generate modulated structures in selected order-parameter directions of isotropy subgroups

## Usage

```python
from pathlib import Path
import phonopy
from phonopy.structure.symmetry import Symmetry
from spgrep_modulation.modulation import Modulation

# Load Phonopy object
path = Path(__file__).resolve().parent.parent / "tests" / "phonopy_mp-2998.yaml.xz"
ph = phonopy.load(path)

# Prepare Modulation class
qpoint = [0.5, 0, 0]  # X point
md = Modulation.with_supercell_and_symmetry_search(
    dynamical_matrix=ph.dynamical_matrix,
    supercell_matrix=[2, 2, 2],
    qpoint=qpoint,
    factor=ph.unit_conversion_factor,
)

# Degenerated imaginary mode
frequency_index = 0
print(f"Frequency (THz): {md.eigvals_to_frequencies(md.eigenspaces[frequency_index][0]):.2f}")
# -> Frequency (THz): -4.88
print(f"Irrep shape: {md.eigenspaces[frequency_index][2].shape}")
# -> Irrep shape: (16, 2, 2)

# Modulated cells corresponding to one-dimensional order-parameter directions of isotropy subgroup
cells = md.get_high_symmetry_modulated_supercells(frequency_index)
for cell in cells:
    symmetry = Symmetry(cell)
    print(f"{symmetry.dataset['international']} (No. {symmetry.dataset['number']})")
# -> Pmma (No. 51) and Cmcm (No. 63)
```

## Installation

```shell
pip install spgrep-modulation
```

```shell
conda create -n modulation python=3.10 pip
conda activate modulation
git clone git@github.com:lan496/spgrep-modulation.git
cd spgrep-modulation
pip install -e .
# pip install -e ".[dev,docs,vis]"
# pre-commit install
```

## License

spgrep-modulation is released under a BSD 3-clause license.

## Development

Document
```shell
sphinx-autobuild docs docs_build
# open localhost:8000 in your browser
```

## Acknowledgements

Some test files `tests/phonopy_mp-*.yaml.xz` are adapted from [phonondb](http://phonondb.mtl.kyoto-u.ac.jp/index.html) under CC BY 4.0.

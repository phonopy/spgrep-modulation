# spgrep-modulation
[![testing](https://github.com/phonopy/spgrep-modulation/actions/workflows/testing.yml/badge.svg)](https://github.com/phonopy/spgrep-modulation/actions/workflows/testing.yml)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Collective atomic modulation analysis with irreducible space-group representation

- Github: <https://github.com/phonopy/spgrep-modulation>
- Document: <https://phonopy.github.io/spgrep-modulation>
- Document (develop): <https://phonopy.github.io/spgrep-modulation/develop/>

## Features

- Calculate representation matrices and irreps formed by phonon eigenmodes
- Calculate isotropy subgroups of irreps of space groups on the fly
- Generate modulated structures in selected order-parameter directions of isotropy subgroups

## Usage

## Installation

```shell
conda create -n spgrep python=3.10 pip
conda activate spgrep
git clone git@github.com:lan496/spgrep-modulation.git
cd spgrep-modulation
pip install -e .
# pip install -e ".[dev,docs]"
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

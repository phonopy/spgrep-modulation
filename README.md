# spgrep-modulation
Collective atomic modulation analysis with irreducible space-group representation

- Document (develop): <https://phonopy.github.io/spgrep-modulation/develop/>

## Installation

```shell
conda create -n spgrep python=3.10 pip
conda activate spgrep
# Install spgrep
git clone git@github.com:spglib/spgrep.git
pushd spgrep
pip install -e .
popd
# Install this package
git clone git@github.com:lan496/spgrep-modulation.git
pushd spgrep-modulation
pip install -e .
# pip install -e ".[dev,docs]"
popd
# pre-commit install
```

## Development

Document
```shell
sphinx-autobuild docs docs_build
# open localhost:8000 in your browser
```

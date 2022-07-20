from pathlib import Path

import pytest

import phonopy
from phonopy import Phonopy


@pytest.fixture(scope="session")
def ph_bto() -> Phonopy:
    path = Path(__file__).resolve().parent / "phonopy_mp-2998.yaml.xz"
    ph = phonopy.load(path)
    return ph

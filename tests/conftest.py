from pathlib import Path

import phonopy
import pytest
from phonopy.api_phonopy import Phonopy


@pytest.fixture(scope="session")
def ph_bto() -> Phonopy:
    path = Path(__file__).resolve().parent / "phonopy_mp-2998.yaml.xz"
    ph = phonopy.load(path)  # type: ignore
    return ph


@pytest.fixture(scope="session")
def ph_mgo() -> Phonopy:
    # Fm-3m (225)
    path = Path(__file__).resolve().parent / "phonopy_mp-1265.yaml.xz"
    ph = phonopy.load(path)  # type: ignore
    return ph


@pytest.fixture(scope="session")
def ph_si_diamond() -> Phonopy:
    # Fd-3m (227)
    path = Path(__file__).resolve().parent / "phonopy_mp-149.yaml.xz"
    ph = phonopy.load(path)  # type: ignore
    return ph


@pytest.fixture(scope="session")
def ph_aln() -> Phonopy:
    # P6_3mc (186)
    # N: (2b, z=0.381)
    # Al: (2b, z=0.498)
    path = Path(__file__).resolve().parent / "phonopy_mp-661.yaml.xz"
    ph = phonopy.load(path)  # type: ignore
    return ph

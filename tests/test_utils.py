import numpy as np
import pytest

from spgrep_modulation.utils import get_commensurate_diagonal_supercell, qr_unique


@pytest.mark.parametrize(
    "qpoint",
    [
        [0, 0, 0],
        [1, 1, 1],
        [1 / 2, 1 / 2, 1 / 2],
        [-1, -1, -1],
        [0.25, 1 / 2, -0.25],
    ],
)
def test_get_commensurate_diagonal_supercell(qpoint):
    """Test qpoints that should fit well into small supercells."""
    supercell = get_commensurate_diagonal_supercell(qpoint)
    for i in supercell @ qpoint:  # should be integers
        rounded = np.rint(i)
        assert np.isclose(rounded - i, 0)


def test_qr_unique():
    a = np.array(
        [
            [1, 0, 1],
            [0, -1, 1],
            [0, 0, -1],
        ]
    )
    q, r = qr_unique(a)
    _, r2 = qr_unique(r)
    assert np.allclose(r, r2)

import numpy as np
import pytest

from spgrep_modulation.utils import (
        get_commensurate_diagonal_supercell,
        qr_unique,
)

@pytest.mark.parametrize(
    "qpoint",
    [
        [0, 0, 0],
        [1, 1, 1],
        [0.33, 0.33, 0.33],
        [1/2, 1/2, 1/2],
        [-1, -1, -1],
        [0.25, 1/2, -0.25],
        [0.59, 0.61, 0.63],
    ],
)
def test_get_commensurate_diagonal_supercell(qpoint):
    supercell = get_commensurate_diagonal_supercell(qpoint)
    for i in supercell @ qpoint: # should be integers
        assert np.isclose(i % 1, 0) or np.isclose(i % 1, 1)    


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

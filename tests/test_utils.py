import numpy as np

from spgrep_modulation.utils import qr_unique


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

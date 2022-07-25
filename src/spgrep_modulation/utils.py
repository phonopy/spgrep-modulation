from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing_extensions import TypeAlias  # for Python<3.10

NDArrayInt: TypeAlias = NDArray[np.intc]
NDArrayFloat: TypeAlias = NDArray[np.float_]
NDArrayComplex: TypeAlias = NDArray[np.complex_]


def qr_unique(a: NDArrayComplex) -> tuple[NDArrayComplex, NDArrayComplex]:
    """Computer QR decomposition
        A = QR,
    where Q is unitary and R is upper triangular and diagonal part of R are chosen to be positive.
    This decomposition is unique if A is full rank.

    Parameters
    ----------
    a: array, (m, n)

    Returns
    -------
    q: array, (m, m)
        unitary
    r: array, (m, n)
        upper triangular and R[i, i] >= 0
    """
    q, r = np.linalg.qr(a)
    for i in range(min(a.shape)):
        if r[i, i] < 0:
            q[:, i] *= -1
            r[i, :] *= -1

    return q, r

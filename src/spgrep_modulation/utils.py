"""Utility functions."""

from __future__ import annotations

from fractions import Fraction
from math import gcd

import numpy as np
from numpy.typing import NDArray
from typing_extensions import TypeAlias  # for Python<3.10

NDArrayInt: TypeAlias = NDArray[np.intc]
NDArrayFloat: TypeAlias = NDArray[np.float64]
NDArrayComplex: TypeAlias = NDArray[np.complex128]


def get_modified_dynamical_matrix(
    dynamical_matrix: NDArrayComplex, scaled_positions: NDArrayFloat, qpoint: NDArrayFloat
):
    r"""Get dynamical matrix phased by only lattice points.

    .. math::
       \\Phi_{\\mu\\mu'}(\\kappa\\kappa'; \\mathbf{q})
            := \\frac{1}{\\sqrt{M_{\\kappa}M_{\\kappa'}}} \\sum_{l'} \\Phi_{\\mu\\mu'}(0\\kappa; l'\\kappa') e^{ i \\mathbf{q} \\cdot \\mathbf{r}(l') }

    Parameters
    ----------
    dynamical_matrix: array, (num_atoms * 3, num_atoms * 3)
        Phonopy's dynamical matrix
    scaled_positions: array, (num_atoms, 3)
    qpoint: array, (3, )

    Returns
    -------
    mdm: array, (num_atoms * 3, num_atoms * 3)
        Dynamical matrix in Maradudin's convention
    """
    num_atoms = len(scaled_positions)
    phase = np.exp(2j * np.pi * np.dot(scaled_positions, qpoint))  # (num_atoms, )
    mdm = (
        dynamical_matrix.reshape(num_atoms, 3, num_atoms, 3)
        * phase[:, None, None, None]
        * np.conj(phase)[None, None, :, None]
    ).reshape(num_atoms * 3, num_atoms * 3)
    return mdm


def qr_unique(a: NDArrayComplex) -> tuple[NDArrayComplex, NDArrayComplex]:
    r"""Compute QR decomposition.

    .. math::
       A = QR,

    where ``Q`` is unitary and ``R`` is upper triangular and diagonal part of ``R`` are chosen to be positive.
    This decomposition is unique if ``A`` is full rank.

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


def get_commensurate_diagonal_supercell(
    qpoint: NDArrayFloat, max_denominator: int = 10, tol: float = 0.01
) -> NDArrayInt:
    """Return minimum diagonal supercell in which ``qpoint`` is commensurate.

    For best accuracy, pass qpoints as fractions, e.g. 1/3 instead of 0.333.

    Parameters
    ----------
    qpoint: point in reciprocal space in reduced coordinates
    max_denominator: largest number of cells in a given direction in real space
    tol: maximum distance from provided qpoint to commensurate qpoint along each reciprocal lattice vector
    """
    assert len(qpoint) == 3
    diag = [0 for _ in range(3)]
    for i, qi in enumerate(qpoint):
        n = Fraction(qi).limit_denominator(max_denominator)
        if not np.isclose(np.abs(n.numerator / n.denominator), np.abs(qi), atol=tol):
            print(
                "Warning: qpoint %s with value %s does not fit with the provided tolerance"
                % (i, qi)
            )
        diag[i] = n.denominator
    return np.diag(diag)


def sample_on_unit_sphere(rng: np.random.Generator, n: int, size: int = 1) -> NDArrayFloat:
    """Return random points from a surface of n-dimensional unit sphere.

    Ref: M. E. Muller, Communications of the ACM 2.4, 19-20 (1959).

    Parameters
    ----------
    rng: numpy's random generator
    n: int
        Number of variables
    size: int, default=1
        Number of points to be sampled

    Returns
    -------
    points: array, (size, n)
    """
    points = rng.standard_normal((size, n))
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    return points


def gcd_on_list(elements: list[int]) -> int:
    """Return greatest common devisor for list of integers."""
    g = elements[0]
    for e in elements[1:]:
        g = gcd(g, e)
    return g


def lcm_on_list(elements: list[int]) -> int:
    """Return least common square for list of integers."""
    lcm = elements[0]
    for e in elements[1:]:
        lcm = lcm * e // gcd(lcm, e)
    return lcm

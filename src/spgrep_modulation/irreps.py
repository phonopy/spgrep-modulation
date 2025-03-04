"""Irreps for phonon."""

from __future__ import annotations

from typing import Literal

import numpy as np
from phonopy.structure.cells import Primitive
from phonopy.structure.symmetry import Symmetry
from spgrep.core import get_spacegroup_irreps_from_primitive_symmetry
from spgrep.representation import project_to_irrep

from spgrep_modulation.utils import NDArrayComplex, NDArrayFloat, NDArrayInt


def project_eigenmode_representation(
    eigenmode_representation: NDArrayComplex,
    primitive: Primitive,
    primitive_symmetry: Symmetry,
    primitive_qpoint: NDArrayFloat,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = 1e-5,
    atol: float = 1e-6,  # Too tight tolerance gives wrong basis vectors...
) -> tuple[list[NDArrayComplex], list[NDArrayComplex], NDArrayInt]:
    """Decompose representation matrices by eigenmodes into irreps.

    Parameters
    ----------
    eigenmode_representation: array, (order, num_atoms, 3, num_atoms, 3)
        Representation matrices formed by phonon eigenmodes
    primitive: phonopy.structure.cells.Primitive
        phonopy's primitive object
    primitive_symmetry: phonopy.structure.symmetry.Symmetry
        phonopy's Symmetry object for primitive cell
    primitive_qpoint: array, (3, )
        q vector in ``primitive``'s dual basis vectors
    method: str, 'Neto' or 'random'
        'Neto': construct irreps from a fixed chain of subgroups of little co-group
        'random': construct irreps by numerically diagonalizing a random matrix commute with regular representation
    rtol: float
        Relative tolerance for comparing float values
    atol: float
        Absolute tolerance to distinguish difference eigenvalues

    Returns
    -------
    all_basis: list
        ``all_basis[i]`` is list of independent basis vectors with (dim, num_atoms, 3) forming ``irreps[i]``.
        Note: phase chosen to be consistent with definition of phonopy's dynamical matrix
    irreps: list of irrep, (little_order, dim, dim)
    mapping_little_group: array[int], (order, )
        list of indices for little group.
    """
    rotations = primitive_symmetry.symmetry_operations["rotations"]
    translations = primitive_symmetry.symmetry_operations["translations"]
    num_atoms = len(primitive)

    small_reps, mapping_little_group = get_spacegroup_irreps_from_primitive_symmetry(
        rotations,
        translations,
        primitive_qpoint,
        method=method,
        rtol=rtol,
        atol=atol,
    )

    modified_basis = []
    irreps = []
    for irrep in small_reps:
        projected = project_to_irrep(
            representation=eigenmode_representation[mapping_little_group].reshape(
                -1, num_atoms * 3, num_atoms * 3
            ),
            irrep=irrep,
            atol=atol,
        )

        if len(projected) == 0:
            continue

        modified_basis.append([basis.reshape(-1, num_atoms, 3) for basis in projected])
        irreps.append(irrep)

    # Apply phase to be consistent with definition of phonopy's dynamical matrix
    phase = np.exp(
        -2j * np.pi * np.dot(primitive.scaled_positions, primitive_qpoint)
    )  # (num_atoms, )
    all_basis = []
    for list_mb in modified_basis:
        basis_for_same_irrep = []
        for mb in list_mb:
            eigenvecs = phase[None, :, None] * mb  # (dim_irrep, num_atoms, 3)
            basis_for_same_irrep.append(eigenvecs)

        all_basis.append(basis_for_same_irrep)

    return all_basis, irreps, mapping_little_group


def get_eigenmode_representation(
    primitive: Primitive,
    primitive_symmetry: Symmetry,
    primitive_qpoint: NDArrayFloat,
) -> NDArrayComplex:
    r"""Compute representation matrix for eigenmodes.

    .. math::
       \Gamma_{\kappa'\mu'; \kappa\mu}^{\mathbf{q}}(g) := \exp \left( -i \mathbf{R}_{g} \mathbf{q} \cdot \mathbf{h}_{g}(\kappa) \right) [\mathbf{R}_{g}]_{\mu'\mu} \delta_{ g\kappa, \kappa' }

    Parameters
    ----------
    primitive: phonopy.structure.cells.Primitive
        phonopy's primitive object
    primitive_symmetry: phonopy.structure.symmetry.Symmetry
        phonopy's Symmetry object for primitive cell
    primitive_qpoint: array, (3, )
        q vector in ``primitive``'s dual basis vectors

    Returns
    -------
    rep: array, (order, num_atoms, 3, num_atoms, 3)
        Representation matrices
    """
    rotations = primitive_symmetry.symmetry_operations["rotations"]
    translations = primitive_symmetry.symmetry_operations["translations"]
    order = len(rotations)
    num_atoms = len(primitive)

    shifts = np.zeros((order, num_atoms, 3))
    for i, (Ri, vi) in enumerate(zip(rotations, translations)):
        perm_i = primitive_symmetry.atomic_permutations[i]
        shifts[i] = (
            primitive.scaled_positions @ Ri.T + vi[None, :] - primitive.scaled_positions[perm_i]
        )

    perm_rep = np.zeros((order, num_atoms, num_atoms), dtype=np.complex128)
    for i, Ri in enumerate(rotations):
        for kappa in range(num_atoms):
            kappa2 = primitive_symmetry.atomic_permutations[i, kappa]
            perm_rep[i, kappa2, kappa] = np.exp(
                -2j * np.pi * np.dot(Ri.T @ primitive_qpoint, shifts[i, kappa])
            )

    # Rotation matrix in cartesian (order, 3, 3)
    A = primitive.cell.T  # column-wise lattice vectors
    Ainv = np.linalg.inv(A)
    rotation_rep = np.array([A @ r @ Ainv for r in rotations], dtype=np.complex128)

    rep = np.einsum("ipq,iab->ipaqb", perm_rep, rotation_rep, optimize="greedy")
    return rep

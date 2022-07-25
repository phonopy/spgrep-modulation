from __future__ import annotations

import numpy as np
from spgrep.group import get_little_group
from spgrep.irreps import enumerate_small_representations
from spgrep.representation import project_to_irrep

from phonopy.harmonic.force_constants import similarity_transformation
from phonopy.structure.cells import Primitive
from phonopy.structure.symmetry import Symmetry
from spgrep_modulation.utils import NDArrayComplex, NDArrayFloat, NDArrayInt


def project_eigenmode_representation(
    eigenmode_representation: NDArrayComplex,
    primitive: Primitive,
    primitive_symmetry: Symmetry,
    primitive_qpoint: NDArrayFloat,
    rtol: float = 1e-5,
    atol: float = 1e-8,
) -> tuple[list[NDArrayComplex], list[NDArrayComplex], NDArrayInt]:
    """
    Parameters
    ----------
    eigenmode_representation: array, (order, num_atoms, 3, num_atoms, 3)
    primitive: Primitive
    primitive_symmetry: Symmetry
    primitive_qpoint: array, (3, )
        q vector in ``primitive``'s dual basis vectors

    Returns
    -------
    basis: list
        basis[i] is list of independent basis vectors with (dim, num_atoms, 3) forming irreps[i].
        Note: phase chosen to be consistent with definition of phonopy's dynamical matrix
    irreps: list of irrep (little_order, dim, dim)
    mapping_little_group:
    """
    rotations = primitive_symmetry.symmetry_operations["rotations"]
    translations = primitive_symmetry.symmetry_operations["translations"]
    num_atoms = len(primitive)

    little_rotations, little_translations, mapping_little_group = get_little_group(
        rotations, translations, primitive_qpoint, atol=atol
    )

    # Compute irreps of little co-group
    little_cogroup_irreps = enumerate_small_representations(
        little_rotations,
        little_translations,
        primitive_qpoint,
        method="Neto",
        rtol=rtol,
        atol=atol,
    )

    modified_basis = []
    irreps = []
    for irrep in little_cogroup_irreps:
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
    basis = []
    for list_mb in modified_basis:
        basis_for_same_irrep = []
        for mb in list_mb:
            eigenvecs = phase[None, :, None] * mb  # (dim_irrep, num_atoms, 3)
            basis_for_same_irrep.append(eigenvecs)

        basis.append(basis_for_same_irrep)

    return basis, irreps, mapping_little_group


def get_eigenmode_representation(
    primitive: Primitive,
    primitive_symmetry: Symmetry,
    primitive_qpoint: NDArrayFloat,
) -> NDArrayComplex:
    """Compute representation matrix for eigenmodes.

    .. math::
       \\Gamma_{\\kappa'\\mu'; \\kappa\\mu}^{\\mathbf{q}}(g) := \\exp \\left( -i \\mathbf{R}_{g} \\mathbf{q} \\cdot \\mathbf{h}_{g}(\\kappa) \\right) [\\mathbf{R}_{g}]_{\\mu'\\mu} \\delta_{ g\\kappa, \\kappa' }

    Parameters
    ----------
    primitive: Primitive
    primitive_symmetry: Symmetry
    primitive_qpoint: array, (3, )
        q vector in ``primitive``'s dual basis vectors

    Returns
    -------
    rep: array, (order, num_atoms, 3, num_atoms, 3)
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

    perm_rep = np.zeros((order, num_atoms, num_atoms), dtype=np.complex_)
    for i, Ri in enumerate(rotations):
        for kappa in range(num_atoms):
            kappa2 = primitive_symmetry.atomic_permutations[i, kappa]
            perm_rep[i, kappa2, kappa] = np.exp(
                -2j * np.pi * np.dot(Ri.T @ primitive_qpoint, shifts[i, kappa])
            )

    # Rotation matrix in cartesian (order, 3, 3)
    rotation_rep = np.array([similarity_transformation(primitive.cell.T, r) for r in rotations])

    rep = np.einsum("ipq,iab->ipaqb", perm_rep, rotation_rep, optimize="greedy")
    return rep

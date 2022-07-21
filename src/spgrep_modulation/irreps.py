import numpy as np

from phonopy.harmonic.force_constants import similarity_transformation
from phonopy.structure.cells import Primitive
from phonopy.structure.symmetry import Symmetry
from spgrep_modulation.utils import NDArrayFloat


def get_eigenmode_representation(
    primitive: Primitive,
    primitive_symmetry: Symmetry,
    primitive_qpoint: NDArrayFloat,
):
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

    rep = np.einsum("ipq,iab->ipaqb", perm_rep, rotation_rep)
    return rep

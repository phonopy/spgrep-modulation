import pytest
from spgrep.group import get_little_group
from spgrep.representation import check_spacegroup_representation, is_unitary

from phonopy.structure.symmetry import Symmetry
from spgrep_modulation.irreps import get_eigenmode_representation


@pytest.mark.parametrize(
    "qpoint,ph_name",
    [
        ([0.0, 0.0, 0.0], "ph_mgo"),  # Gamma point
        ([0.0, 0.5, 0.5], "ph_mgo"),  # M point
        ([0.0, 0.0, 0.0], "ph_bto"),  # Gamma point
        ([0.0, 0.5, 0.5], "ph_bto"),  # M point
    ],
)
def test_eigenmode_representation(request, qpoint, ph_name):
    ph = request.getfixturevalue(ph_name)

    primitive = ph.dynamical_matrix.primitive
    primitive_symmetry = Symmetry(cell=primitive)

    rep = get_eigenmode_representation(primitive, primitive_symmetry, qpoint)

    little_rotations, little_translations, mapping = get_little_group(
        primitive_symmetry.symmetry_operations["rotations"],
        primitive_symmetry.symmetry_operations["translations"],
        qpoint,
    )

    num_atoms = len(primitive)
    # Check if `rep` is unitary representation
    assert is_unitary(rep.reshape(-1, num_atoms * 3, num_atoms * 3))
    # Check if `rep` preserves multiplication for little group
    assert check_spacegroup_representation(
        little_rotations,
        little_translations,
        qpoint,
        rep[mapping].reshape(-1, num_atoms * 3, num_atoms * 3),
    )

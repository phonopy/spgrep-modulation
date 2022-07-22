import pytest
from spgrep.group import get_little_group
from spgrep.irreps import enumerate_small_representations
from spgrep.representation import check_spacegroup_representation, is_unitary

from phonopy.structure.symmetry import Symmetry
from spgrep_modulation.irreps import (
    get_eigenmode_representation,
    project_eigenmode_representation,
)


@pytest.mark.parametrize(
    "ph_name,qpoint,num_irreps,num_basis",
    [
        ("ph_si_diamond", [0, 0, 0], 2, [1, 1]),
        ("ph_si_diamond", [0.5, 0, 0.5], 3, [1, 1, 1]),  # X point in primitive
        ("ph_mgo", [0, 0, 0], 1, [2]),
        ("ph_mgo", [0.5, 0, 0.5], 2, [2, 2]),  # X point in primitive
    ],
)
def test_project_eigenmode_representation(request, ph_name, qpoint, num_irreps, num_basis):
    ph = request.getfixturevalue(ph_name)

    primitive = ph.dynamical_matrix.primitive
    primitive_symmetry = Symmetry(cell=primitive)

    rep = get_eigenmode_representation(primitive, primitive_symmetry, qpoint)
    basis, irreps, mapping = project_eigenmode_representation(
        rep, primitive, primitive_symmetry, qpoint
    )
    assert len(irreps) == num_irreps
    assert [len(b) for b in basis] == num_basis


@pytest.mark.parametrize(
    "ph_name,qpoint",
    [
        ("ph_mgo", [0.0, 0.0, 0.0]),  # Gamma point
        ("ph_mgo", [0.0, 0.5, 0.5]),  # M point
        ("ph_bto", [0.0, 0.0, 0.0]),  # Gamma point
        ("ph_bto", [0.0, 0.5, 0.5]),  # M point
        ("ph_si_diamond", [0, 0, 0]),
        ("ph_si_diamond", [0.5, 0, 0.5]),  # X point in primitive
    ],
)
def test_eigenmode_representation(request, ph_name, qpoint):
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

    irreps = enumerate_small_representations(little_rotations, little_translations, qpoint)
    for irrep in irreps:
        assert check_spacegroup_representation(
            little_rotations,
            little_translations,
            qpoint,
            irrep,
        )

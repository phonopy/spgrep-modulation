import numpy as np
import pytest
from phonopy.structure.symmetry import Symmetry
from spgrep.group import get_cayley_table
from spgrep.pointgroup import pg_dataset
from spgrep.utils import is_integer_array

from spgrep_modulation.isotropy import (
    IsotropyEnumerator,
    enumerate_point_subgroup,
    enumerate_point_subgroup_naive,
    get_translational_subgroup,
    search_compliment,
)
from spgrep_modulation.modulation import Modulation


@pytest.mark.parametrize(
    "qpoint",
    [
        ([0, 0, 0]),
        ([0.5, 0, 0]),
        ([0, 0.5, 0]),
        ([0.5, 0.5, 0]),
        ([0.5, 0.5, 0.5]),
        ([1 / 3, 0, 0]),
        ([2 / 3, 0, 0]),
        ([1 / 3, 1 / 3, 0]),
        ([2 / 3, 2 / 3, 0]),
        ([1 / 3, 2 / 3, 0]),
        ([1 / 3, 1 / 3, 1 / 3]),
        ([1 / 3, 1 / 3, 2 / 3]),
        ([1 / 3, 2 / 3, 2 / 3]),
        ([2 / 3, 2 / 3, 2 / 3]),
        ([0.5, 1 / 3, 2 / 3]),
        ([0, 0, 0.5]),
    ],
)
def test_get_translational_subgroup(qpoint):
    transformation = get_translational_subgroup(qpoint)
    assert np.linalg.det(transformation) > 0
    assert is_integer_array(transformation @ qpoint)


def test_enumerate_point_subgroup():
    pointgroup = pg_dataset["4/mmm"][0]
    table = get_cayley_table(np.array(pointgroup))
    flags = [True for _ in range(len(pointgroup))]
    subgroups_actual = enumerate_point_subgroup(table, flags, return_conjugacy_class=False)
    subgroups_expect = enumerate_point_subgroup_naive(table, flags)
    assert len(subgroups_actual) == len(subgroups_expect)


def test_compliments():
    rotations = np.array(
        [
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
            [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
            [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
            [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
            [[0, -1, 0], [-1, 0, 0], [0, 0, 1]],
            [[0, -1, 0], [-1, 0, 0], [0, 0, -1]],
            [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
        ]
    )
    translations = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.5],
        ]
    )

    transformation = np.diag([1, 1, 2])
    table = get_cayley_table(rotations)
    point_subgroup_indices = enumerate_point_subgroup(table, [True for _ in rotations])
    search_compliment(rotations, translations, point_subgroup_indices[-1], transformation, table)
    for indices in point_subgroup_indices:
        result = search_compliment(rotations, translations, indices, transformation, table)
        if indices == [0, 7]:
            assert not result


# Ref: Table 1 of "Isotropy Subgroups of the 230 Crystallographic Space Groups" (p.349)
@pytest.mark.parametrize(
    "ph_name,qpoint,dimension,freq_idx,numbers_expect",
    [
        # Pm-3m (No. 221), two-dimensional irrep: X_5+
        # Size=2 isotropy subgroups for X_5+
        #   One-dimensional OPD: No. 51, 63
        #   Two-dimensional OPD: No. 11
        ("ph_bto", [0, 0.5, 0], [2, 2, 2], 0, [11, 51, 63]),
        # P6_3mc (186), two-dimensional irrep: K_3
        ("ph_aln", [1 / 3, 1 / 3, 0.0], [3, 3, 2], 0, None),
    ],
)
def test_isotropy_subgroup(request, ph_name, qpoint, dimension, freq_idx, numbers_expect):
    ph = request.getfixturevalue(ph_name)

    md = Modulation.with_supercell_and_symmetry_search(
        dynamical_matrix=ph.dynamical_matrix,
        supercell_matrix=dimension,
        qpoint=qpoint,
        factor=ph.unit_conversion_factor,
    )

    _, _, irrep = md.eigenspaces[freq_idx]
    assert irrep.shape[1] == 2

    ie = IsotropyEnumerator(
        md.little_rotations,
        md.little_translations,
        qpoint,
        irrep,
    )

    numbers_actual = []
    maximal_displacement = 0.11
    for opd in ie.order_parameter_directions:
        # For more than two-dimensional OPD, choose the first vector
        amplitudes = np.abs(opd)[0]
        arguments = np.angle(opd)[0]
        modulation = md.get_modulated_supercell_and_modulation(
            freq_idx, amplitudes, arguments, return_cell=False
        )
        scaled_modulation = maximal_displacement / np.max(np.abs(modulation)) * modulation
        cell = md.apply_modulation_to_supercell(scaled_modulation)
        symmetry = Symmetry(cell)

        numbers_actual.append(symmetry.dataset["number"])

    if numbers_expect:
        assert set(numbers_actual) == set(numbers_expect)

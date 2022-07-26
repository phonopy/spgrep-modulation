import numpy as np
import pytest

from phonopy import Phonopy
from phonopy.phonon.modulation import Modulation as PhonopyModulation
from spgrep_modulation.modulation import Modulation


@pytest.mark.parametrize(
    "ph_name,qpoint,dimension",
    [
        ("ph_bto", [0, 0, 0], [1, 1, 1]),  # Gamma
        ("ph_bto", [0, 0, 0.5], [2, 2, 2]),  # X point
        ("ph_bto", [0, 0.5, 0.5], [2, 2, 2]),  # M point
        ("ph_mgo", [0, 0, 0], [1, 1, 1]),  # Gamma
        ("ph_mgo", [0.5, 0, 0.5], [2, 2, 2]),  # X point in primitive
        ("ph_si_diamond", [0, 0, 0], [1, 1, 1]),  # Gamma
        ("ph_si_diamond", [0.5, 0, 0.5], [2, 2, 2]),  # X point in primitive
        ("ph_aln", [0.0, 0.0, 0.0], [1, 1, 1]),
        ("ph_aln", [0.0, 0.0, 0.5], [2, 2, 2]),  # A point
        ("ph_aln", [1 / 2, 0.0, 0.0], [2, 2, 2]),  # M point
        ("ph_aln", [1 / 3, 1 / 3, 0.0], [3, 3, 2]),  # K point
    ],
)
def test_symmetry_adapted_eigenmodes(request, ph_name, qpoint, dimension):
    ph = request.getfixturevalue(ph_name)

    md = Modulation.with_supercell_and_symmetry_search(
        dynamical_matrix=ph.dynamical_matrix,
        supercell_matrix=dimension,
        qpoint=qpoint,
        factor=ph.unit_conversion_factor,
    )

    # Check if each mode is truly eigenvector of dynamical matrix
    dm = md.dynamical_matrix.dynamical_matrix
    num_atoms = len(md.primitive)
    for eigval, modes in md.eigenspaces:
        actual = np.einsum("ij,kj->ki", dm, modes.reshape(-1, num_atoms * 3), optimize="greedy")
        expect = eigval * modes.reshape(-1, num_atoms * 3)
        assert np.allclose(actual, expect, atol=1e-5)


def test_regression(ph_bto: Phonopy):
    # Supercell dimension with respect to the primitive cell, dtype='intc', shape=(3, ), (3, 3), (9, )
    dimension = [2, 2, 2]

    qpoint = [0.0, 0.5, 0.5]  # M point, non-degenerated
    band_index = 0
    amplitude = 2.0
    argument = 0.0
    phonon_modes = [
        [
            qpoint,  # q-point in reduced coordinates
            band_index,
            amplitude * np.sqrt(len(ph_bto.primitive)) / 2,
            np.degrees(argument),
        ],
    ]

    phmd = PhonopyModulation(
        dynamical_matrix=ph_bto.dynamical_matrix,
        dimension=dimension,
        phonon_modes=phonon_modes,
        factor=ph_bto.unit_conversion_factor,
    )
    phmd.run()
    modulation_phonopy = phmd._u[0]

    md = Modulation.with_supercell_and_symmetry_search(
        dynamical_matrix=ph_bto.dynamical_matrix,
        supercell_matrix=dimension,
        qpoint=qpoint,
        factor=ph_bto.unit_conversion_factor,
    )
    modulated_cell, modulation = md.get_modulated_supercell_and_modulation(
        frequency_index=band_index,
        amplitudes=[amplitude],
        arguments=[argument],
    )

    # Check modulations up to U(1)
    indices = np.nonzero(modulation_phonopy)
    phase_diff = (0.5 * modulation[indices[0][0], indices[1][0]]) / modulation_phonopy[
        indices[0][0], indices[1][0]
    ]
    assert np.allclose(0.5 * modulation, modulation_phonopy * phase_diff)
    assert np.isclose(np.abs(phase_diff), 1.0)

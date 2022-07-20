import numpy as np

from phonopy import Phonopy
from phonopy.phonon.modulation import Modulation as PhonopyModulation
from phonopy.structure.atoms import PhonopyAtoms
from spgrep_modulation.modulation import Modulation


def compare_cells_with_order(cell: PhonopyAtoms, cell_ref: PhonopyAtoms, symprec=1e-5):
    """Compare two cells with the same orders of positions."""
    np.testing.assert_allclose(cell.cell, cell_ref.cell, atol=symprec)
    compare_positions_with_order(cell.scaled_positions, cell_ref.scaled_positions, cell.cell)
    np.testing.assert_array_equal(cell.numbers, cell_ref.numbers)
    np.testing.assert_allclose(cell.masses, cell_ref.masses, atol=symprec)
    if cell.magnetic_moments is None:
        assert cell_ref.magnetic_moments is None
    else:
        np.testing.assert_allclose(cell.magnetic_moments, cell_ref.magnetic_moments, atol=symprec)


def compare_positions_with_order(pos, pos_ref, lattice, symprec=1e-5):
    """Compare two lists of positions and orders."""
    diff = pos - pos_ref
    diff -= np.rint(diff)
    dist = (np.dot(diff, lattice) ** 2).sum(axis=1)
    assert (dist < symprec).all()


def test_regression(ph_bto: Phonopy):
    # Supercell dimension with respect to the primitive cell, dtype='intc', shape=(3, ), (3, 3), (9, )
    dimension = [2, 2, 2]

    # qpoint = [0., 0., 0.5]  # X point, degenerated as dim=2
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
    modulated_cell_phonopy = phmd.get_modulated_supercells()[0]

    md = Modulation.with_supercell_and_symmetry_search(
        dynamical_matrix=ph_bto.dynamical_matrix,
        supercell_matrix=dimension,
        qpoint=qpoint,
        factor=ph_bto.unit_conversion_factor,
    )
    modulated_cell = md.get_modulated_supercell(
        frequency_index=band_index,
        amplitudes=[amplitude],
        arguments=[argument + np.pi],  # shift phase to match with phonopy's modulation
    )

    compare_cells_with_order(modulated_cell, modulated_cell_phonopy)

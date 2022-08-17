import pytest

from spgrep_modulation.isotropy import IsotropyEnumerator
from spgrep_modulation.modulation import Modulation


@pytest.mark.parametrize(
    "ph_name,qpoint,dimension,freq_idx",
    [
        ("ph_bto", [0, 0, 0.5], [2, 2, 2], 0),  # X point
        # ("ph_aln", [1 / 3, 1 / 3, 0.0], [3, 3, 2], 0),  # K point
    ],
)
def test_isotropy_subgroup(request, ph_name, qpoint, dimension, freq_idx):
    ph = request.getfixturevalue(ph_name)

    md = Modulation.with_supercell_and_symmetry_search(
        dynamical_matrix=ph.dynamical_matrix,
        supercell_matrix=dimension,
        qpoint=qpoint,
        factor=ph.unit_conversion_factor,
    )

    # Two-dimensional irrep
    _, _, irrep = md.eigenspaces[freq_idx]
    assert irrep.shape[1] == 2

    ie = IsotropyEnumerator(
        md.little_rotations,
        md.little_translations,
        qpoint,
        irrep,
    )
    print(ie.order_parameter_directions)

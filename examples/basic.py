from pathlib import Path

import phonopy
from phonopy.structure.symmetry import Symmetry

from spgrep_modulation.modulation import Modulation

# Load Phonopy object
path = Path(__file__).resolve().parent.parent / "tests" / "phonopy_mp-2998.yaml.xz"
ph = phonopy.load(path)  # type: ignore

# Prepare Modulation class
qpoint = [0.5, 0, 0]  # X point
md = Modulation.with_supercell_and_symmetry_search(
    dynamical_matrix=ph.dynamical_matrix,
    supercell_matrix=[2, 2, 2],
    qpoint=qpoint,
    factor=ph.unit_conversion_factor,
)

# Degenerated imaginary mode
frequency_index = 0
print(f"Frequency (THz): {md.eigvals_to_frequencies(md.eigenspaces[frequency_index][0]):.2f}")
# -> Frequency (THz): -4.88
print(f"Irrep shape: {md.eigenspaces[frequency_index][2].shape}")
# -> Irrep shape: (16, 2, 2)

# Modulated cells corresponding to one-dimensional order-parameter directions of isotropy subgroup
cells = md.get_high_symmetry_modulated_supercells(frequency_index)
for cell in cells:
    symmetry = Symmetry(cell)
    print(f"{symmetry.dataset['international']} (No. {symmetry.dataset['number']})")
# -> Pmma (No. 51) and Cmcm (No. 63)

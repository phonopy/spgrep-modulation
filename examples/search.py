from __future__ import annotations

import json
from collections import deque
from dataclasses import dataclass
from logging import INFO, Formatter, StreamHandler, basicConfig, getLogger
from pathlib import Path

import click
import numpy as np
import seekpath
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT
from ase.constraints import UnitCellFilter
from ase.optimize import BFGS
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.structure.symmetry import Symmetry
from phonopy.units import VaspToTHz
from pymatgen.analysis.prototypes import AflowPrototypeMatcher
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.phonopy import get_pmg_structure

from spgrep_modulation.modulation import Modulation
from spgrep_modulation.utils import NDArrayFloat, get_commensurate_diagonal_supercell


@dataclass
class TrialCell:
    cell: PhonopyAtoms
    relaxed: PhonopyAtoms | None = None
    relaxed_energy_per_atom: float | None = None
    parent: PhonopyAtoms | None = None
    qpoint: NDArrayFloat | None = None
    frequency_index: int | None = None


class ModulationSearch:
    def __init__(
        self,
        calc: Calculator,
        supercell_matrix: list[list[int]] = [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
        maximal_displacement: float = 0.11,
        max_size: int = 512,
        symmetry_tolerance: float = 1e-2,  # Rough symprec to idealize modulated cells
    ) -> None:
        self._calc = calc
        self._supercell_matrix = supercell_matrix
        self._maximal_displacement = maximal_displacement
        self._max_size = max_size
        self._symmetry_tolerance = symmetry_tolerance
        self._structure_matcher = StructureMatcher(
            ltol=1e-7, stol=1e-7, scale=False
        )  # Tighter tolerance than default

        # Logger
        self._logger = getLogger(__name__)
        basicConfig(level=INFO)
        log_fmt = Formatter(
            "%(asctime)s %(name)s %(lineno)d [%(levelname)s][%(funcName)s] %(message)s "
        )
        handler = StreamHandler()
        handler.setLevel("INFO")
        handler.setFormatter(log_fmt)

    @property
    def calc(self) -> Calculator:
        return self._calc

    @property
    def supercell_matrix(self) -> list[list[int]]:
        return self._supercell_matrix

    @property
    def maximal_displacement(self) -> float:
        return self._maximal_displacement

    @property
    def max_size(self) -> int:
        return self._max_size

    @property
    def symmetry_tolerance(self) -> float:
        return self._symmetry_tolerance

    def run_recursive(
        self, initial: TrialCell, max_depth: int = 4
    ) -> tuple[list[TrialCell], list[int]]:
        found = [
            initial,
        ]
        found_pmg = [
            get_pmg_structure(initial.cell),
        ]
        parents = [
            -1,
        ]

        que = deque()  # type: ignore
        que.append((initial, 0, 0))
        while len(que) > 0:
            parent, idx, depth = que.pop()
            if depth > max_depth:
                continue

            self._logger.info(f"parent_id={idx}, depth={depth + 1}")
            children = self.run(parent)

            for child in children:
                child_pmg = get_pmg_structure(child.cell)
                if any([self._structure_matcher.fit(child_pmg, other) for other in found_pmg]):
                    continue

                found.append(child)
                found_pmg.append(child_pmg)
                parents.append(idx)
                que.append((child, len(found), depth + 1))

        return found, parents

    def run(self, initial: TrialCell) -> list[TrialCell]:
        # Relax cell
        relaxed_cell, energy = self._get_relaxed_cell_and_energy(
            initial.cell,
            mask=[True, True, True, True, True, True],
        )
        initial.relaxed = relaxed_cell
        initial.relaxed_energy_per_atom = energy / len(relaxed_cell)

        # Calculate harmonic force constants
        ph = self._get_phonon(relaxed_cell, self.supercell_matrix)
        # Create modulation generators and search instable modes
        modulators_imag_freqs = self._search_instable_modes(relaxed_cell, ph)

        next_trialcells = []
        for md, imag_freqs in modulators_imag_freqs:
            for idx in imag_freqs:
                next_cells_idx = md.get_modulated_supercells(
                    idx, self.maximal_displacement, self.max_size
                )
                selected = self._pick_and_refine_high_symmetry_cells(next_cells_idx)

                for cell in selected:
                    next_trial = TrialCell(
                        cell=cell,
                        parent=relaxed_cell,
                        qpoint=md.qpoint,
                        frequency_index=idx,
                    )
                    next_trialcells.append(next_trial)

        if len(next_trialcells) == 0:
            self._logger.info("End this branch.")

        return next_trialcells

    def _get_energy_and_forces(self, cell: PhonopyAtoms) -> tuple[float, NDArrayFloat]:
        atoms = Atoms(
            symbols=cell.symbols, scaled_positions=cell.scaled_positions, cell=cell.cell, pbc=True
        )
        atoms.set_calculator(self.calc)
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        return energy, forces

    def _get_relaxed_cell_and_energy(
        self, cell: PhonopyAtoms, mask: list[bool]
    ) -> tuple[PhonopyAtoms, float]:
        atoms = Atoms(
            symbols=cell.symbols, scaled_positions=cell.scaled_positions, cell=cell.cell, pbc=True
        )
        atoms.set_calculator(self.calc)

        ucf = UnitCellFilter(atoms, mask=mask)
        dyn = BFGS(ucf)
        dyn.run()
        energy = atoms.get_potential_energy()

        relaxed_cell = PhonopyAtoms(
            symbols=atoms.symbols,
            scaled_positions=atoms.get_scaled_positions(),
            cell=atoms.cell,
        )
        self._logger.info(f"Finish to relax structure: energy={energy/len(atoms):.4f} eV/atom")
        return relaxed_cell, energy

    def _get_phonon(self, cell: PhonopyAtoms, supercell_matrix, distance: float = 0.03) -> Phonopy:
        ph = Phonopy(
            cell,
            supercell_matrix=supercell_matrix,
            primitive_matrix="auto",
            factor=VaspToTHz,
        )
        ph.generate_displacements(distance=distance)

        # Prepare force sets
        force_sets = []
        for supercell in ph.supercells_with_displacements:
            _, forces = self._get_energy_and_forces(supercell)
            force_sets.append(forces)
        ph.forces = force_sets
        ph.produce_force_constants()

        self._logger.info("Finish to calculate force constants.")
        return ph

    def _search_instable_modes(
        self, cell: PhonopyAtoms, ph: Phonopy
    ) -> list[tuple[Modulation, list[int]]]:
        band_path = seekpath.get_path(cell.totuple())
        point_coords = band_path["point_coords"]

        modulators_imag_freqs = []
        for label, qpoint in point_coords.items():
            dimension = get_commensurate_diagonal_supercell(qpoint)
            # Skip general points
            if np.isclose(np.linalg.det(dimension), 0):
                continue

            md = Modulation.with_supercell_and_symmetry_search(
                dynamical_matrix=ph.dynamical_matrix,
                supercell_matrix=dimension,
                qpoint=qpoint,
                factor=ph.unit_conversion_factor,
            )

            imag_freqs = []
            for idx, (eigval, eigvecs, _) in enumerate(md.eigenspaces):
                freq = md.eigvals_to_frequencies(eigval)
                # Search imaginary mode
                if freq < -md.degeneracy_tolerance:
                    imag_freqs.append(idx)
                    degeneracy = len(eigvecs)
                    self._logger.info(
                        f"Instable mode at {label}={qpoint} with degeneracy={degeneracy}."
                    )

            if len(imag_freqs) == 0:
                continue

            modulators_imag_freqs.append((md, imag_freqs))

        return modulators_imag_freqs

    def _pick_and_refine_high_symmetry_cells(
        self, cells: list[PhonopyAtoms]
    ) -> list[PhonopyAtoms]:
        # Ref: https://github.com/atztogo/cogue/blob/master/cogue/phonon/modulation.py
        max_num_ops = 0
        selected = []
        selected_pmg = []
        symbols = []
        for cell in cells:
            symmetry = Symmetry(cell, symprec=self.symmetry_tolerance)
            refined = PhonopyAtoms(
                scaled_positions=symmetry.dataset["std_positions"],
                cell=symmetry.dataset["std_lattice"],
                numbers=symmetry.dataset["std_types"],
            )
            num_ops = len(symmetry.symmetry_operations["rotations"])

            if num_ops > max_num_ops:
                max_num_ops = num_ops
                selected = [
                    refined,
                ]
                selected_pmg = [
                    get_pmg_structure(refined),
                ]
                symbols = [
                    symmetry.get_international_table(),
                ]
            elif num_ops == max_num_ops:
                pmg_structure = get_pmg_structure(refined)
                if any(
                    [self._structure_matcher.fit(pmg_structure, other) for other in selected_pmg]
                ):
                    continue
                selected.append(refined)
                selected_pmg.append(pmg_structure)
                symbols.append(symmetry.get_international_table())

        self._logger.info(
            f"Pick {len(selected)} high symmetry cells (num_ops={max_num_ops}) out of {len(cells)}."
        )
        self._logger.info(f"Space-group types of picked cells: {symbols}")
        return selected


@click.command()
@click.option(
    "--output", required=False, default="debug", type=str, help="Root directory for saved files"
)
def main(output):
    calc = EMT()  # Effective medium potential for Al, Ni, Cu, Pd, Ag, Pt and Au
    supercell_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]

    ms = ModulationSearch(
        calc,
        supercell_matrix,
        max_size=256,
        maximal_displacement=0.15,
        symmetry_tolerance=1e-3,
    )

    # Simple cubic structure
    a = 2.4
    initial_cell = PhonopyAtoms(
        symbols=["Cu"],
        scaled_positions=[[0.0, 0.0, 0.0]],
        cell=[[a, 0, 0], [0, a, 0], [0, 0, a]],
    )

    found, parents = ms.run_recursive(TrialCell(cell=initial_cell))

    root = Path(output)
    root.mkdir(exist_ok=True)

    apm = AflowPrototypeMatcher()

    prototypes = []
    for i, tcell in enumerate(found):
        structure = get_pmg_structure(tcell.cell)

        prototype = apm.get_prototypes(structure)
        prototype_tag = ""
        if prototype:
            prototype_tag = prototype[0]["tags"]["aflow"]
        prototypes.append(prototype_tag)
        print(f"Energy: {tcell.relaxed_energy_per_atom:.4f} eV/atom, Prototype: {prototype_tag}")

        path = root / f"structure_{i}.cif"
        structure.to(fmt="cif", filename=path, symprec=1e-5)

    result = {
        "energies": [tcell.relaxed_energy_per_atom for tcell in found],
        "paths": [(src, dst) for dst, src in enumerate(parents)],
        "prototypes": prototypes,
    }
    with open(root / "result.json", "w") as f:
        json.dump(result, f, indent=2)


if __name__ == "__main__":
    main()

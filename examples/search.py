from __future__ import annotations

import json
from abc import ABC, abstractmethod
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
from ase.filters import FrechetCellFilter
from ase.optimize import LBFGSLineSearch
from phonopy.api_phonopy import Phonopy
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
    idx: int
    parent_idx: int
    cell: PhonopyAtoms
    relaxed: PhonopyAtoms | None = None
    relaxed_energy_per_atom: float | None = None
    qpoint: NDArrayFloat | None = None
    frequency_index: int | None = None


class AbstractModulationSearch(ABC):
    def __init__(
        self,
        calc: Calculator,
        supercell_matrix: list[list[int]] = [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
        maximal_displacement: float = 0.11,
        max_size: int = 512,
        symmetry_tolerance: float = 1e-2,  # Rough symprec to idealize modulated cells
        ltol: float = 1e-8,
        stol: float = 1e-8,
        angle_tol: float = -1,
    ) -> None:
        self._calc = calc
        self._supercell_matrix = supercell_matrix
        self._maximal_displacement = maximal_displacement
        self._max_size = max_size
        self._symmetry_tolerance = symmetry_tolerance
        self._structure_matcher = StructureMatcher(
            ltol=ltol,
            stol=stol,
            angle_tol=angle_tol,
            scale=False,
        )  # Tighter tolerance than default

        # Global variable to track index of cells
        self._max_idx = 1

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

    def run_recursive(self, initial_cell: PhonopyAtoms, max_depth: int = 4) -> list[TrialCell]:
        initial = TrialCell(
            idx=0,
            cell=initial_cell,
            parent_idx=-1,
        )
        self._max_idx = 1

        found = []
        found_pmg = []

        que: deque[tuple[TrialCell, int]] = deque()
        que.append((initial, 0))
        while len(que) > 0:
            parent, depth = que.pop()
            if depth > max_depth:
                self._logger.info(f"Stop child of parent={parent.idx} due to max depth limit")
                continue

            self._logger.info(f"parent={parent.idx}, depth={depth}")

            # Show symmetry
            symmetry = Symmetry(parent.cell, symprec=1e-5)
            self._logger.info(
                f"{symmetry.dataset['international']} ({symmetry.dataset['number']})"
            )

            # Relax cell
            relaxed_cell, energy = self._get_relaxed_cell_and_energy(
                parent.cell,
                mask=[True, True, True, True, True, True],
            )
            parent.relaxed = relaxed_cell
            parent.relaxed_energy_per_atom = energy / len(relaxed_cell)

            # Now, entries of `parent` are filled
            found.append(parent)
            parent_pmg = get_pmg_structure(parent.cell)
            found_pmg.append(parent_pmg)

            children = self.run(parent)

            for child in children:
                child_pmg = get_pmg_structure(child.cell)
                if any([self._structure_matcher.fit(child_pmg, other) for other in found_pmg]):
                    self._logger.info("Detect duplicates")
                    continue

                que.appendleft((child, depth + 1))

        return found

    @abstractmethod
    def run(self, initial: TrialCell) -> list[TrialCell]:
        pass

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

        filter = FrechetCellFilter(atoms, mask=mask)
        dyn = LBFGSLineSearch(filter)
        dyn.run()
        energy = atoms.get_potential_energy()

        relaxed_cell = PhonopyAtoms(
            symbols=atoms.symbols,
            scaled_positions=atoms.get_scaled_positions(),
            cell=atoms.cell,
        )
        self._logger.info(f"Finish to relax structure: energy={energy / len(atoms):.4f} eV/atom")
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


class BaseModulationSearch(AbstractModulationSearch):
    def run(self, initial: TrialCell) -> list[TrialCell]:
        # Calculate harmonic force constants
        ph = self._get_phonon(initial.relaxed, self.supercell_matrix)
        # Create modulation generators and search instable modes
        modulators_imag_freqs = self._search_instable_modes(initial.relaxed, ph)

        next_trialcells = []
        for md, imag_freqs in modulators_imag_freqs:
            for idx in imag_freqs:
                next_cells_idx = md.get_modulated_supercells(
                    idx, self.maximal_displacement, self.max_size
                )
                selected = self._pick_and_refine_high_symmetry_cells(next_cells_idx)

                for cell in selected:
                    next_trial = TrialCell(
                        idx=self._max_idx,
                        parent_idx=initial.idx,
                        cell=cell,
                        qpoint=md.qpoint,
                        frequency_index=idx,
                    )
                    next_trialcells.append(next_trial)
                    self._max_idx += 1

        if len(next_trialcells) == 0:
            self._logger.info("End this branch.")

        return next_trialcells

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


class IsotropyModulationSearch(AbstractModulationSearch):
    def run(self, initial: TrialCell) -> list[TrialCell]:
        # Calculate harmonic force constants
        ph = self._get_phonon(initial.relaxed, self.supercell_matrix)
        # Create modulation generators and search instable modes
        modulators_imag_freqs = self._search_instable_modes(initial.relaxed, ph)

        next_trialcells = []
        for md, imag_freqs in modulators_imag_freqs:
            for idx in imag_freqs:
                next_cells_idx = md.get_high_symmetry_modulated_supercells(
                    idx, self.maximal_displacement
                )
                self._logger.info(f"{len(next_cells_idx)} high-symmetry modulations are detected.")
                for cell in next_cells_idx:
                    next_trial = TrialCell(
                        idx=self._max_idx,
                        parent_idx=initial.idx,
                        cell=cell,
                        qpoint=md.qpoint,
                        frequency_index=idx,
                    )
                    next_trialcells.append(next_trial)
                    self._max_idx += 1

        if len(next_trialcells) == 0:
            self._logger.info("End this branch.")

        return next_trialcells


@click.command()
@click.option(
    "--output", required=False, default="debug", type=str, help="Root directory for saved files"
)
@click.option(
    "--show_prototype", required=False, default=False, help="Attach prototype information if true"
)
def main(output, show_prototype):
    calc = EMT()  # Effective medium potential for Al, Ni, Cu, Pd, Ag, Pt and Au
    supercell_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]

    ms = IsotropyModulationSearch(
        calc,
        supercell_matrix,
        ltol=1e-6,
        stol=1e-6,
        angle_tol=2.0,
    )

    # Simple cubic structure
    a = 2.4
    initial_cell = PhonopyAtoms(
        symbols=["Cu"],
        scaled_positions=[[0.0, 0.0, 0.0]],
        cell=[[a, 0, 0], [0, a, 0], [0, 0, a]],
    )

    found = ms.run_recursive(initial_cell, max_depth=2)

    root = Path(output)
    root.mkdir(exist_ok=True)

    apm = AflowPrototypeMatcher()

    prototypes = []
    for i, tcell in enumerate(found):
        structure = get_pmg_structure(tcell.cell)

        if show_prototype:
            prototype = apm.get_prototypes(structure)
            prototype_tag = ""
            if prototype:
                prototype_tag = prototype[0]["tags"]["aflow"]
            prototypes.append(prototype_tag)
            print(
                f"Energy: {tcell.relaxed_energy_per_atom:.4f} eV/atom, Prototype: {prototype_tag}"
            )
        else:
            print(f"Energy: {tcell.relaxed_energy_per_atom:.4f} eV/atom")

        path = root / f"structure_{i}.cif"
        structure.to(fmt="cif", filename=path, symprec=1e-5)

    result = {
        "energies": [tcell.relaxed_energy_per_atom for tcell in found],
        "paths": [(tcell.parent_idx, tcell.idx) for tcell in found],
    }
    if show_prototype:
        result["prototypes"] = prototypes

    with open(root / "result.json", "w") as f:
        json.dump(result, f, indent=2)


if __name__ == "__main__":
    main()

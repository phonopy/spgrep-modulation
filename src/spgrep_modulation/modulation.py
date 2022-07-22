from __future__ import annotations

from warnings import warn

import numpy as np

from phonopy.harmonic.dynamical_matrix import DynamicalMatrix, DynamicalMatrixNAC
from phonopy.phonon.degeneracy import degenerate_sets, get_eigenvectors
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.structure.cells import (
    Primitive,
    Supercell,
    get_supercell,
    is_primitive_cell,
    shape_supercell_matrix,
)
from phonopy.structure.symmetry import Symmetry
from phonopy.units import VaspToTHz
from spgrep_modulation.utils import NDArrayComplex, NDArrayFloat, NDArrayInt


class Modulation:
    """
    Parameters
    ----------
    primitive_symmetry:
    supercell:
    dynamical_matrix:
    qpoint: array, (3, )
    nac_q_direction:
    factor:
    degeneracy_tolerance: float
        Absolute tolerance to groupby phonon frequencies in ``factor`` unit
    """

    def __init__(
        self,
        primitive_symmetry: Symmetry,
        supercell: Supercell,
        dynamical_matrix: DynamicalMatrix | DynamicalMatrixNAC,
        qpoint: NDArrayFloat,
        nac_q_direction: NDArrayFloat | None = None,
        factor: float = VaspToTHz,
        degeneracy_tolerance: float = 1e-5,
    ) -> None:
        # Check to be commensurate
        if not np.allclose(np.remainder(supercell.supercell_matrix.T @ qpoint, 1), 0):
            warn(f"Given qpoint={qpoint} is not commensurate with supercell.")
        self._qpoint = qpoint

        self._primitive_symmetry = primitive_symmetry
        if not is_primitive_cell(primitive_symmetry.symmetry_operations["rotations"]):
            raise RuntimeError("Set primitive cell.")

        self._supercell = supercell
        self._dynamical_matrix = dynamical_matrix
        self._nac_q_direction = nac_q_direction
        self._factor = factor

        self._supercell_size = np.abs(np.around(np.linalg.det(self._supercell.supercell_matrix)))

        # Diagonalize dynamical matrix if not performed yet
        eigvals, eigvecs = get_eigenvectors(
            qpoint,
            self._dynamical_matrix,
            ddm=None,  # Not used
            perturbation=None,
            derivative_order=None,
            nac_q_direction=self._nac_q_direction,
        )

        # Group eigenvecs by frequencies
        self._eigenspaces = []
        self._frequencies = self._eigvals_to_frequencies(eigvals)
        deg_sets = degenerate_sets(self._frequencies, cutoff=degeneracy_tolerance)
        for indices in deg_sets:
            self._eigenspaces.append((eigvals[indices[0]], [eigvecs[:, idx] for idx in indices]))

        # TODO: Construct irreps with `qpoint` here

    @property
    def primitive_symmetry(self) -> Symmetry:
        return self._primitive_symmetry

    @property
    def primitive(self) -> Primitive:
        return self._dynamical_matrix.primitive

    @property
    def supercell(self) -> Supercell:
        return self._supercell

    @property
    def unit_conversion_factor(self) -> float:
        return self._factor

    @property
    def frequencies(self) -> NDArrayFloat:
        return self._frequencies

    def get_modulated_supercell_and_modulation(
        self,
        frequency_index: int,
        amplitudes: list[float],
        arguments: list[float],
    ) -> tuple[PhonopyAtoms, NDArrayComplex]:
        # Adapted from phonopy
        _, eigvecs = self._eigenspaces[frequency_index]

        # Generate modulation
        modulation = np.zeros((len(self._supercell), 3), dtype=np.complex_)
        for eigvec, amplitude, argument in zip(eigvecs, amplitudes, arguments):
            modulation += self._get_displacements(eigvec.reshape(-1, 3), amplitude, argument)

        # Apply modulation to supercell
        lattice = self._supercell.cell
        positions = self._supercell.positions
        positions += np.real(modulation) / 2
        scaled_positions = np.dot(positions, np.linalg.inv(lattice))
        scaled_positions = np.remainder(scaled_positions, 1)
        cell = self._supercell.copy()
        cell.scaled_positions = scaled_positions

        return cell, modulation

    @classmethod
    def with_supercell_and_symmetry_search(
        cls,
        dynamical_matrix: DynamicalMatrix | DynamicalMatrixNAC,
        supercell_matrix: NDArrayInt,
        qpoint: NDArrayFloat,
        nac_q_direction: NDArrayFloat | None = None,
        factor: float = VaspToTHz,
        symprec: float = 1e-5,
    ) -> Modulation:
        primitive = dynamical_matrix.primitive
        primitive_symmetry = Symmetry(cell=primitive, symprec=symprec)

        supercell_matrix_3x3 = shape_supercell_matrix(supercell_matrix)
        supercell = get_supercell(
            unitcell=dynamical_matrix.primitive,
            supercell_matrix=supercell_matrix_3x3,
            is_old_style=False,  # Use Smith normal form
        )

        return cls(
            primitive_symmetry=primitive_symmetry,
            supercell=supercell,
            dynamical_matrix=dynamical_matrix,
            qpoint=qpoint,
            nac_q_direction=nac_q_direction,
            factor=factor,
        )

    def _eigvals_to_frequencies(self, eigvals: NDArrayComplex) -> NDArrayFloat:
        # Adapted from phonopy
        e = np.array(eigvals).real
        return np.sqrt(np.abs(e)) * np.sign(e) * self._factor

    def _get_displacements(
        self,
        reshaped_eigvec: NDArrayComplex,  # (primitive, 3)
        amplitude: float,
        argument: float,  # in radian
    ) -> NDArrayComplex:
        # Adapted from phonopy
        m = self._supercell.masses
        s2uu_map = [self._supercell.u2u_map[x] for x in self._supercell.s2u_map]
        spos = self._supercell.scaled_positions
        dim = self._supercell.supercell_matrix
        coeffs = np.exp(2j * np.pi * np.dot(np.dot(spos, dim.T), self._qpoint)) / np.sqrt(m)
        # Change denominator from phonopy's modulation: len(m) -> self._supercell_size
        # to adapt equation in the formulation page.
        u = reshaped_eigvec[s2uu_map, :] * coeffs[:, None] / np.sqrt(self._supercell_size)

        # Omit calibration of phase in phonopy's modulation for now
        u *= amplitude * np.exp(1j * argument)

        return u

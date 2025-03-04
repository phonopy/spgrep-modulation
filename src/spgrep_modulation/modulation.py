"""Modulation class."""

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

from spgrep_modulation.irreps import (
    get_eigenmode_representation,
    project_eigenmode_representation,
)
from spgrep_modulation.isotropy import IsotropyEnumerator
from spgrep_modulation.utils import (
    NDArrayComplex,
    NDArrayFloat,
    NDArrayInt,
    get_modified_dynamical_matrix,
    qr_unique,
    sample_on_unit_sphere,
)


class Modulation:
    """Generate modulated cells based on dynamical matrix of phonon."""

    def __init__(
        self,
        primitive_symmetry: Symmetry,
        supercell: Supercell,
        dynamical_matrix: DynamicalMatrix | DynamicalMatrixNAC,
        qpoint: NDArrayFloat,
        nac_q_direction: NDArrayFloat | None = None,
        factor: float = VaspToTHz,
        degeneracy_tolerance: float = 1e-4,
        seed: int = 0,
    ) -> None:
        """Generate modulated cells based on dynamical matrix of phonon.

        Parameters
        ----------
        primitive_symmetry: phonopy.structure.symmetry.Symmetry
            phonopy's Symmetry object for primitive cell
        supercell: phonopy.structure.cells.Supercell
            phonopy's supercell object to apply modulation
        dynamical_matrix: phonopy.harmonic.dynamical_matrix.DynamicalMatrix
            phonopy's dynamical matrix object
        qpoint: array, (3, )
        nac_q_direction:
            Used to calculate dynamical matrix around Gamma point if specified
        factor: float
            Unit conversion
        degeneracy_tolerance: float
            Absolute tolerance to groupby phonon frequencies in ``factor`` unit
        seed: int
            Seed value for sampling modulation from degenerated modes
        """
        # Check to be commensurate
        if not np.allclose(np.remainder(supercell.supercell_matrix.T @ qpoint, 1), 0):
            warn(f"Given qpoint={qpoint} is not commensurate with supercell.")
        self._qpoint = qpoint

        self._primitive_symmetry = primitive_symmetry
        rotations = primitive_symmetry.symmetry_operations["rotations"]
        translations = primitive_symmetry.symmetry_operations["translations"]
        if not is_primitive_cell(rotations):
            raise RuntimeError("Set primitive cell.")

        self._supercell = supercell
        self._dynamical_matrix = dynamical_matrix
        self._nac_q_direction = nac_q_direction
        self._factor = factor
        self._degeneracy_tolerance = degeneracy_tolerance
        self._rng = np.random.default_rng(seed=seed)

        self._supercell_size = np.abs(np.around(np.linalg.det(self._supercell.supercell_matrix)))

        # Diagonalize dynamical matrix if not performed yet
        eigvals, _ = get_eigenvectors(
            qpoint,
            self._dynamical_matrix,
            ddm=None,  # Not used
            perturbation=None,
            derivative_order=None,
            nac_q_direction=self._nac_q_direction,
        )

        # Group eigenvecs by frequencies
        self._eigenspaces, mapping_little_group = self._group_eigenspaces(eigvals)

        # Little group
        self._little_rotations = rotations[mapping_little_group]
        self._little_translations = translations[mapping_little_group]

    def _group_eigenspaces(self, eigvals):
        # Construct irreps with `qpoint`
        rep = get_eigenmode_representation(self.primitive, self.primitive_symmetry, self.qpoint)
        all_basis, irreps, mapping_little_group = project_eigenmode_representation(
            rep,
            self.primitive,
            self.primitive_symmetry,
            self.qpoint,
            atol=self.degeneracy_tolerance,
        )

        # Modified dynamical matrix
        mdm = get_modified_dynamical_matrix(
            self.dynamical_matrix.dynamical_matrix, self.primitive.scaled_positions, self.qpoint
        )

        eigenspaces = []
        num_atoms = len(self.primitive)
        phase = np.exp(
            2j * np.pi * np.dot(self.primitive.scaled_positions, self.qpoint)
        )  # (num_atoms, )
        for list_basis, irrep in zip(all_basis, irreps):
            # Block-diagonalize modified dynamical matrix
            list_modified_basis = [basis * phase[None, :, None] for basis in list_basis]
            F_irrep = np.concatenate(
                [mb.reshape(-1, num_atoms * 3) for mb in list_modified_basis]
            ).T  # (num_atoms * 3, len(list_basis) * dim_irrep)
            mdm_irrep = (
                np.conj(F_irrep.T) @ mdm @ F_irrep
            )  # (len(list_basis) * dim_irrep, len(list_basis) * dim_irrep)
            if not np.allclose(mdm_irrep, np.conj(mdm_irrep).T):
                warn("Block-diagonalized modified dynamical matrix is not Hermitian.")

            eigvals_irrep, eigvecs_irrep = np.linalg.eigh(mdm_irrep)

            # Group by eigenvalues
            frequencies_irrep = self.eigvals_to_frequencies(eigvals_irrep)
            cutoff = self.degeneracy_tolerance
            max_groupby_degeneracy_trials = 8
            for _ in range(max_groupby_degeneracy_trials):
                # TODO: binary search for choosing cutoff value
                deg_sets_irrep = degenerate_sets(frequencies_irrep, cutoff=cutoff)
                if len(deg_sets_irrep) < len(list_basis):
                    cutoff /= 10  # Tighten tolerance
                elif len(deg_sets_irrep) == len(list_basis):
                    break

            if len(deg_sets_irrep) != len(list_basis):
                warn(
                    f"If no accidental degeneracy happens, number of unique eigenvalues corresponding to an irrep ({len(list_basis)}) should be equal to degeneracy of the irrep ({len(deg_sets_irrep)})."
                )

            dim_irrep = irrep.shape[1]
            for indices in deg_sets_irrep:
                # If no accidental degeneracy happens, deg_eigval == dim_irrep, len(deg_sets_irrep) == len(list_basis)
                deg_eigval = len(indices)
                if (not np.any(np.isclose(eigvals, eigvals_irrep[indices[0]]))) or (
                    deg_eigval != dim_irrep
                ):
                    warn(f"Inconsistent eigenvalue: {eigvals_irrep[indices[0]]}, {eigvals}")

                space = np.array(
                    [eigvecs_irrep[:, idx] for idx in indices]
                )  # (deg_eigval=dim_irrep, len(list_basis) * dim_irrep)
                # QR decomposition of column-wise vectors gives Gram-Schmidt orthonormalized vectors in column wise.
                space = qr_unique(space.T)[0].T
                # Go back eigenvectors with phonopy's convention
                space_phonopy = np.einsum(
                    "k,kmp,ip->ikm",
                    np.conj(phase),
                    F_irrep.reshape(num_atoms, 3, len(list_basis) * dim_irrep),
                    space,
                    optimize="greedy",
                )  # (deg_eigval=dim_irrep, num_atoms, 3)

                # Representation matrix
                irrep_eigmodes = np.einsum(
                    "nsq,kqp,msp->knm",
                    np.conj(space).reshape(deg_eigval, len(list_basis), dim_irrep),
                    irrep,  # (little_order, dim_irrep, dim_irrep)
                    space.reshape(deg_eigval, len(list_basis), dim_irrep),
                    optimize="greedy",
                )  # (little_order, dim_irrep, dim_irrep)

                eigenspaces.append((eigvals_irrep[indices[0]], space_phonopy, irrep_eigmodes))

        # Sort by eigenvalue for regression test
        sorted_eigenspaces = sorted(eigenspaces, key=lambda e: e[0])

        return sorted_eigenspaces, mapping_little_group

    @property
    def dynamical_matrix(self) -> DynamicalMatrix | DynamicalMatrixNAC:
        """Return Phonopy's dynamical matrix."""
        return self._dynamical_matrix

    @property
    def primitive_symmetry(self) -> Symmetry:
        """Return Phonopy's symmetry for primitive cell."""
        return self._primitive_symmetry

    @property
    def little_rotations(self) -> NDArrayInt:
        """Return rotations of little group at ``qpoint``."""
        return self._little_rotations

    @property
    def little_translations(self) -> NDArrayFloat:
        """Return translations of little group at ``qpoint``."""
        return self._little_translations

    @property
    def primitive(self) -> Primitive:
        """Return Phonopy's primitive cell."""
        return self._dynamical_matrix.primitive

    @property
    def supercell(self) -> Supercell:
        """Return Phonopy's supercell without modulations."""
        return self._supercell

    @property
    def qpoint(self) -> NDArrayFloat:
        """Return qpoint."""
        return self._qpoint

    @property
    def unit_conversion_factor(self) -> float:
        """Return unit conversion for frequencies."""
        return self._factor

    @property
    def degeneracy_tolerance(self) -> float:
        """Return absolute tolerance to detect degenerated bands."""
        return self._degeneracy_tolerance

    @property
    def eigenspaces(self) -> list[tuple[float, NDArrayComplex, NDArrayComplex]]:
        r"""Return list of ``(eigenvalue, degenerated eigenvectors, irrep formed by the eigenvectors)`` of dynamical matrix.

        ``frequency_index`` is the index of ``eigenspaces``.
        Each ``frequency_index`` corresponds to a eigenvalue of the dynamical matrix, :math:`\omega_{\alpha s}^{2}`.
        Here, :math:`\alpha` specifies irrep up to unitary equivalence and :math:`s` distinguishes the irrep among equivalent irreps.
        ``eigenspaces[frequency_index]`` is a tuple with

        * ``eigenvalue``: float, eigenvalue of the dynamical matrix :math:`\omega_{\alpha s}^{2}`
        * ``eigenvectors``: array, (dim, num_atoms, 3), ``eigenvectors[lambda, kappa, :]`` is an eigenvector :math:`\mathbf{e}(\kappa; \mathbf{q}\alpha s \lambda)`.
        * ``irrep``: array, (little_order, dim, dim), representation matrices of the irrep :math:`\tilde{\Gamma}^{\mathbf{q}\alpha s}`

        """
        return self._eigenspaces

    def get_modulated_supercell_and_modulation(
        self,
        frequency_index: int,
        amplitudes: list[float],
        arguments: list[float],
        return_cell: bool = True,
    ) -> tuple[PhonopyAtoms, NDArrayComplex] | NDArrayComplex:
        """Return modulated cell and modulation with given amplitudes and arguments.

        Parameters
        ----------
        frequency_index: int
            Index of considered eigenmodes in ``Modulation.eigenspaces``.
            See :func:`Modulation.eigenspaces`.
        amplitudes: list[float], (dim, )
        arguments: list[float], (dim, )
            in radian
        return_cell: bool = True
        """
        # Adapted from phonopy
        _, eigvecs, _ = self.eigenspaces[frequency_index]

        # Generate modulation
        modulation = np.zeros((len(self.supercell), 3), dtype=np.complex128)
        for eigvec, amplitude, argument in zip(eigvecs, amplitudes, arguments):
            modulation += self._get_displacements(eigvec.reshape(-1, 3), amplitude, argument)

        if return_cell:
            cell = self.apply_modulation_to_supercell(modulation)
            return cell, modulation
        else:
            return modulation

    def get_modulated_supercells(
        self,
        frequency_index: int,
        maximal_displacement: float = 0.11,
        max_size: int = 1,
    ) -> list[PhonopyAtoms]:
        """Return modulated cells with randomly sampled arguments.

        Parameters
        ----------
        frequency_index: int
            Index of considered eigenmodes in ``Modulation.eigenspaces``.
            See :func:`Modulation.eigenspaces`.
        maximal_displacement: int
            Amplitude of modulation is chosen so that the maximal displacement in returned modulation is equal to ``maximal_displacement``.
            Same unit as ``Modulation.primitive.cell``.
        max_size: int
            Number of returned modulated structures when degenerated eigenmodes are specified.

        Returns
        -------
        cells: list of PhonopyAtoms
        """
        _, eigvecs, _ = self.eigenspaces[frequency_index]
        degeneracy = len(eigvecs)

        modulations = []
        if degeneracy == 1:
            modulation = self.get_modulated_supercell_and_modulation(
                frequency_index, amplitudes=[1.0], arguments=[0.0], return_cell=False
            )
            modulations.append(modulation)
        elif degeneracy > 1:
            points = np.zeros((max_size, 2 * degeneracy))
            # One of order parameters can be chosen as pure real because of unitary arbitrariness of eigenvectors
            points[:, :-1] = sample_on_unit_sphere(self._rng, 2 * degeneracy - 1, size=max_size)
            order_params = points.reshape(max_size, degeneracy, 2)
            order_params = order_params[:, :, 0] + order_params[:, :, 1] * 1.0j

            amplitudes = np.abs(order_params)
            arguments = np.angle(order_params)
            for amp, arg in zip(amplitudes, arguments):
                modulation = self.get_modulated_supercell_and_modulation(
                    frequency_index, amp, arg, return_cell=False
                )
                modulations.append(modulation)

        cells = []
        for modulation in modulations:
            # Scale modulation to so that its maximal displacement is equal to ``maximal_displacement``
            scaled_modulation = maximal_displacement / np.max(np.abs(modulation)) * modulation
            cell = self.apply_modulation_to_supercell(scaled_modulation)
            cells.append(cell)

        return cells

    def get_high_symmetry_modulated_supercells(
        self,
        frequency_index: int,
        maximal_displacement: float = 0.11,
    ) -> list[PhonopyAtoms]:
        """Return modulated cells corresponding to one-dimensional order-parameter directions of an isotropy subgroup.

        Parameters
        ----------
        frequency_index: int
            Index of considered eigenmodes in ``Modulation.eigenspaces``.
            See :func:`Modulation.eigenspaces`.
        maximal_displacement: int
            Amplitude of modulation is chosen so that the maximal displacement in returned modulation is equal to ``maximal_displacement``.
            Same unit as ``Modulation.primitive.cell``.

        Returns
        -------
        cells: list of PhonopyAtoms
        """
        _, _, irrep = self.eigenspaces[frequency_index]
        ie = IsotropyEnumerator(
            self.little_rotations,
            self.little_translations,
            self.qpoint,
            irrep,
        )

        cells = []

        # Choose isotropy subgroups with one-dimensional order-parameter directions
        for opd in ie.order_parameter_directions:
            if len(opd) > 1:
                continue
            amplitudes = np.abs(opd)[0]
            arguments = np.angle(opd)[0]
            modulation = self.get_modulated_supercell_and_modulation(
                frequency_index, amplitudes, arguments, return_cell=False
            )
            scaled_modulation = maximal_displacement / np.max(np.abs(modulation)) * modulation
            cell = self.apply_modulation_to_supercell(scaled_modulation)
            cells.append(cell)

        return cells

    @classmethod
    def with_supercell_and_symmetry_search(
        cls,
        dynamical_matrix: DynamicalMatrix | DynamicalMatrixNAC,
        supercell_matrix: NDArrayInt,
        qpoint: NDArrayFloat,
        nac_q_direction: NDArrayFloat | None = None,
        factor: float = VaspToTHz,
        degeneracy_tolerance: float = 1e-4,
        symprec: float = 1e-5,
        seed: int = 0,
    ) -> Modulation:
        """Return Modulation class with supercell matrix."""
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
            degeneracy_tolerance=degeneracy_tolerance,
            factor=factor,
            seed=seed,
        )

    def eigvals_to_frequencies(self, eigvals: NDArrayComplex) -> NDArrayFloat:
        """Convert eigenvalue to frequency."""
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

    def apply_modulation_to_supercell(self, modulation: NDArrayComplex) -> PhonopyAtoms:
        """Apply modulation to supercell."""
        lattice = self.supercell.cell
        positions = self.supercell.positions
        positions += np.real(modulation) / 2
        scaled_positions = np.dot(positions, np.linalg.inv(lattice))
        scaled_positions = np.remainder(scaled_positions, 1)
        cell = self.supercell.copy()
        cell.scaled_positions = scaled_positions
        return cell

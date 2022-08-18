from __future__ import annotations

from fractions import Fraction
from queue import Queue

import numpy as np
from hsnf import row_style_hermite_normal_form
from hsnf.integer_system import solve_modular_integer_linear_system
from spgrep.group import get_cayley_table, get_identity_index, get_inverse_index
from spgrep.representation import project_to_irrep
from spgrep.utils import contain_space, is_integer_array

from spgrep_modulation.utils import (
    NDArrayComplex,
    NDArrayFloat,
    NDArrayInt,
    lcm_on_list,
)


class IsotropyEnumerator:
    """
    maximal_isotropy_subgroups: list of indices
    order_parameter_directions:
    """

    def __init__(
        self,
        little_rotations: NDArrayInt,
        little_translations: NDArrayFloat,
        qpoint: NDArrayFloat,
        small_rep: NDArrayComplex,
        max_denominator: int = 100,
        atol: float = 1e-6,
    ) -> None:
        self._little_rotations = np.array(little_rotations)
        self._little_translations = np.array(little_translations)
        self._qpoint = np.array(qpoint)
        self._small_rep = small_rep
        self._max_denominator = max_denominator
        self._atol = atol

        self._maximal_isotropy_subgroups, self._order_parameter_directions = self._initialize()

    @property
    def little_rotations(self) -> NDArrayInt:
        return self._little_rotations

    @property
    def little_translations(self) -> NDArrayFloat:
        return self._little_translations

    @property
    def qpoint(self) -> NDArrayFloat:
        return self._qpoint

    @property
    def small_rep(self) -> NDArrayComplex:
        return self._small_rep

    @property
    def maximal_isotropy_subgroups(self) -> list[list[int]]:
        return self._maximal_isotropy_subgroups

    @property
    def order_parameter_directions(self) -> list[NDArrayComplex]:
        return self._order_parameter_directions

    def _initialize(self):
        transformation = self._get_translational_subgroup()

        # Point group operations that preserve translational subgroup of isotropy subgroup
        preserve_sublattice = [False for _ in range(len(self.little_rotations))]
        for i, R in enumerate(self.little_rotations):
            RM = np.linalg.inv(transformation) @ R @ transformation
            if is_integer_array(RM, atol=self._atol):
                preserve_sublattice[i] = True

        table = get_cayley_table(self.little_rotations)
        point_subgroup_indices = self._enumerate_point_subgroup(preserve_sublattice, table)
        isotropy_subgroups = []
        order_parameter_directions = []
        for indices in point_subgroup_indices:
            if not self._is_space_subgroup(indices, transformation, table):
                continue

            directions = self._determine_order_parameter_directions(indices)
            if len(directions) == 0:
                continue

            # Group by subspace spanned by order-parameter directions
            found = False
            for i, other in enumerate(order_parameter_directions):
                if (len(other) == len(directions)) and contain_space(
                    other, directions, atol=self._atol
                ):
                    isotropy_subgroups[i].append(indices)
                    found = True
                    break
            if not found:
                order_parameter_directions.append(directions)
                isotropy_subgroups.append([indices])

        # Choose largest isotropy with the same order-parameter directions (subduction criterion)
        maximal_isotropy_subgroups = []
        for subgroups in isotropy_subgroups:
            idx = np.argmax([len(subgroup) for subgroup in subgroups])
            maximal_isotropy_subgroups.append(subgroups[idx])

        return maximal_isotropy_subgroups, order_parameter_directions

    def _get_translational_subgroup(self):
        # Basis vectors of a sublattice formed by translation that preserve order parameter
        elements = [
            Fraction(qi).limit_denominator(self._max_denominator).denominator for qi in self.qpoint
        ]
        lcm = lcm_on_list(elements)
        A = np.around(self.qpoint * lcm).astype(int)[None, :]
        # basis vectors for k.t = 0
        basis_orth, _ = solve_modular_integer_linear_system(A, np.zeros(1), lcm)
        # basis vector parallel to k
        basis_para = np.around(self.qpoint * lcm).astype(int) * lcm

        # transformation @ mathbb{Z}^{3} forms sublattice
        transformation = np.vstack([basis_orth, basis_para])
        if np.linalg.det(transformation) < 0:
            transformation[0, :] *= -1
        transformation, _ = row_style_hermite_normal_form(transformation)

        return transformation

    def _enumerate_point_subgroup_naive(self, preserve_sublattice: list[bool]):
        """Not used, but for checking number with the bit-DP version."""
        table = get_cayley_table(self.little_rotations)
        order = len(table)
        ret = []
        for bits in range(1, 1 << order):
            elements = self._decode_bits(bits, order)
            if not all([preserve_sublattice[idx] for idx in elements]):
                continue

            if self._is_subgroup(elements, table):
                ret.append(bits)

        return ret

    def _enumerate_point_subgroup(
        self,
        preserve_sublattice: list[bool],
        table: NDArrayInt,
    ) -> list[list[int]]:
        order = len(table)
        identity = get_identity_index(table)
        # Represent choice of elements by bit array
        st = {1 << identity}
        for i in range(order):
            if (not preserve_sublattice[i]) or (i == identity):
                continue
            if (1 << i) in st:
                # Already visited
                continue

            next_st = set()
            for bits in st:
                elements = self._decode_bits(bits, order)
                assert self._is_subgroup(elements, table)
                generated = self._traverse(elements + [i], identity, table)
                next_st.add(sum(1 << idx for idx in set(generated)))

            st = st.union(next_st)

        # Group by conjugacy classes
        found = set()
        ret = []
        for bits in sorted(st):
            if bits in found:
                continue
            elements = self._decode_bits(bits, order)
            ret.append(elements)
            for i in range(order):
                if not preserve_sublattice[i]:
                    continue
                inv = get_inverse_index(table, i)
                conj = [table[inv, table[idx, i]] for idx in elements]
                found.add(sum(1 << idx for idx in set(conj)))

        assert len(found) == len(st)
        return ret

    def _is_space_subgroup(
        self,
        indices: list[int],
        transformation: NDArrayInt,
        table: NDArrayInt,
    ) -> bool:
        """Check translational parts"""
        # TODO: check additional translation for non-symmorphic space groups
        Minv = np.linalg.inv(transformation)
        for i in indices:
            Ri = self.little_rotations[i]
            vi = self.little_translations[i]
            inv_i = get_inverse_index(table, i)
            inv_Ri = np.linalg.inv(Ri)
            for j in indices:
                vj = self._little_translations[j]
                k = table[inv_i, j]
                disp = vj - inv_Ri @ vi - self.little_translations[k]
                if not is_integer_array(Minv @ disp, atol=self._atol):
                    return False

        return True

    def _determine_order_parameter_directions(
        self, isotropy_subgroup: list[int]
    ) -> list[NDArrayComplex]:
        """Determine order-parameter directions of given isotropy subgroup.

        Returns
        -------
        directions: (num_directions, small_rep_dim)
        """
        subduced_rep = self.small_rep[isotropy_subgroup]
        id_rep = np.ones((len(isotropy_subgroup), 1, 1))
        directions = project_to_irrep(subduced_rep, id_rep, atol=self._atol)
        return np.array(directions).reshape(-1, self.small_rep.shape[1])

    @staticmethod
    def _decode_bits(bits: int, order: int) -> list[int]:
        elements = [idx for idx in range(order) if (bits >> idx) & 1 == 1]
        return elements

    @staticmethod
    def _is_subgroup(elements: list[int], table: NDArrayInt) -> bool:
        subtable = table[elements][:, elements]
        for i in range(len(subtable)):
            if (set(subtable[i]) != set(elements)) or (set(subtable[:, i]) != set(elements)):
                return False
        return True

    @staticmethod
    def _traverse(
        generators: list[int],
        identity: int,
        table: NDArrayInt,
    ) -> list[int]:
        subgroup = set()
        que = Queue()  # type: ignore
        que.put(identity)

        while not que.empty():
            g = que.get()
            if g in subgroup:
                continue
            subgroup.add(g)

            for h in generators:
                que.put(table[g, h])

        return sorted(list(subgroup))

from __future__ import annotations

from fractions import Fraction
from itertools import product
from queue import Queue

import numpy as np
from hsnf import row_style_hermite_normal_form, smith_normal_form
from hsnf.integer_system import solve_integer_linear_system
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
        transformation = get_translational_subgroup(self.qpoint, self._max_denominator)
        assert is_integer_array(transformation @ self.qpoint, atol=self._atol)

        # Point group operations that preserve translational subgroup of isotropy subgroup
        preserve_sublattice = [False for _ in range(len(self.little_rotations))]
        for i, R in enumerate(self.little_rotations):
            RM = np.linalg.inv(transformation) @ R @ transformation
            if is_integer_array(RM, atol=self._atol):
                preserve_sublattice[i] = True

        table = get_cayley_table(self.little_rotations)
        point_subgroup_indices = enumerate_point_subgroup(table, preserve_sublattice)
        isotropy_subgroups = []
        order_parameter_directions = []
        for indices in point_subgroup_indices:
            if not search_compliment(
                self.little_rotations, self.little_translations, indices, transformation, table
            ):
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


def get_translational_subgroup(qpoint: NDArrayFloat, max_denominator: int = 100):
    """Return transformation matrix of the following sublattice:
    Let `t` be a lattice point of the sublattice. Then, np.dot(t, qpoint) is integer.
    """
    if isinstance(qpoint, list):
        qpoint = np.array(qpoint)

    # Basis vectors of a sublattice formed by translation that preserve order parameter
    elements = [Fraction(qi).limit_denominator(max_denominator).denominator for qi in qpoint]
    lcm = lcm_on_list(elements)
    A = np.around(qpoint * lcm).astype(int)[None, :]

    all_basis = []
    # Three vectors [lcm, 0, 0], [0, lcm, 0], [0, 0, lcm] are trivially general solutions.
    all_basis.append(lcm * np.eye(3, dtype=int))

    # Solve `A @ t = lcm`
    # Since GCD of A is 1, this integer linear system always has a special solution.
    basis_and_special = solve_integer_linear_system(A, np.array([lcm]))
    if basis_and_special:
        all_basis.append(basis_and_special[0])
        all_basis.append(basis_and_special[1])

    # transformation @ mathbb{Z}^{3} forms sublattice
    transformation, _ = row_style_hermite_normal_form(np.vstack(all_basis))
    # Take only three basis vectors
    transformation = transformation[:3]

    if np.linalg.det(transformation) < 0:
        transformation[0, :] *= -1

    return transformation


def enumerate_point_subgroup(
    table: NDArrayInt, preserve_sublattice: list[bool], return_conjugacy_class: bool = True
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
            elements = _decode_bits(bits, order)
            assert _is_subgroup(elements, table)
            generated = traverse(elements + [i], identity, table)
            next_st.add(sum(1 << idx for idx in set(generated)))

        st = st.union(next_st)

    if not return_conjugacy_class:
        subgroups = []
        for bits in sorted(st):
            subgroups.append(_decode_bits(bits, order))
        return subgroups

    # Group by conjugacy classes
    found = set()
    ret = []
    for bits in sorted(st):
        if bits in found:
            continue
        elements = _decode_bits(bits, order)
        ret.append(elements)
        for i in range(order):
            if not preserve_sublattice[i]:
                continue
            inv = get_inverse_index(table, i)
            conj = [table[inv, table[idx, i]] for idx in elements]
            found.add(sum(1 << idx for idx in set(conj)))

    assert len(found) == len(st)
    return ret


def enumerate_point_subgroup_naive(table, preserve_sublattice: list[bool]):
    order = len(table)
    ret = []
    for bits in range(1, 1 << order):
        elements = _decode_bits(bits, order)
        if not all([preserve_sublattice[idx] for idx in elements]):
            continue

        if _is_subgroup(elements, table):
            ret.append(bits)

    return ret


def _decode_bits(bits: int, order: int) -> list[int]:
    elements = [idx for idx in range(order) if (bits >> idx) & 1 == 1]
    return elements


def _is_subgroup(elements: list[int], table: NDArrayInt) -> bool:
    subtable = table[elements][:, elements]
    for i in range(len(subtable)):
        if (set(subtable[i]) != set(elements)) or (set(subtable[:, i]) != set(elements)):
            return False
    return True


def traverse(
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


def search_compliment(
    rotations,
    translations,
    indices: list[int],
    transformation: NDArrayInt,
    table: NDArrayInt,
    atol: float = 1e-6,
) -> bool:
    """Return true if some adjusted translations enable given point group and translational subgroup correspond to a space group."""
    # Search generators
    identity = get_identity_index(table)
    generators = _search_generators(indices, identity, table)

    # Lattice points in sublattice formed by `transformation`
    # See https://lan496.github.io/dsenum/supercell.html
    snf, L, R = smith_normal_form(transformation)  # snf = L @ transformation @ R
    invariant_factors = tuple(snf.diagonal())
    points = []
    for factor in product(*[range(f) for f in invariant_factors]):
        point = np.linalg.inv(L) @ np.array(factor)
        points.append(np.around(point).astype(int))

    Minv = np.linalg.inv(transformation)
    for trials in product(points, repeat=len(generators)):
        adjusted_translations = [None for _ in translations]
        adjusted_translations[identity] = np.array([0, 0, 0])
        for g, trial in zip(generators, trials):
            adjusted_translations[g] = translations[g] + trial

        # Determine adjusted translations during traversal
        subgroup: set = set()  # type: ignore
        que = Queue()  # type: ignore
        que.put((identity, adjusted_translations[identity]))
        while not que.empty():
            g, trns = que.get()
            if g in subgroup:
                continue
            subgroup.add(g)
            if adjusted_translations[g] is None:
                adjusted_translations[g] = trns

            for h in generators:
                trns_gh = rotations[g] @ adjusted_translations[h] + adjusted_translations[g]
                que.put((table[g, h], trns_gh))

        # Check consistency
        ok = True
        for g in indices:
            if not is_integer_array(adjusted_translations[g] - translations[g]):
                ok = False
                break
        if not ok:
            continue

        for g, h in product(indices, repeat=2):
            gh = table[g, h]
            additional = (
                rotations[g] @ adjusted_translations[h]
                + adjusted_translations[g]
                - adjusted_translations[gh]
            )
            if not is_integer_array(Minv @ additional, atol=atol):
                ok = False
                break
        if ok:
            return True

    return False


def _search_generators(indices: list[int], identity: int, table: NDArrayInt) -> list[int]:
    generators = []
    for i in indices:
        if i == identity:
            continue
        generators.append(i)
        subgroup = traverse(generators, identity, table)
        if len(subgroup) == len(indices):
            break
    return generators

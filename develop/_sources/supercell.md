# Supercell construction

## Klassengleiche subgroup

Let {math}`L(\mathcal{G})` be a lattice of space group {math}`\mathcal{G}`,
```{math}
    L(\mathcal{G}) = \left\{ \mathbf{t} \mid (\mathbf{I}, \mathbf{t}) \in \mathcal{T}(\mathcal{G}) \right\}.
```

Let {math}`\mathcal{P} := \mathcal{P}(\mathcal{G})` be a point group of space group {math}`\mathcal{G}`.
Lattice {math}`L (\subset L(\mathcal{G}))` is called {math}`\mathcal{P}`-invariant if {math}`\mathcal{P}` acts faithfully on {math}`L` with the following action, {math}`\mathcal{P} \times L \ni (\mathbf{R}, \mathbf{t}) \mapsto \mathbf{Rt} \in L`.

Theorem 5.2 in Ref. {cite}`Eick:js0047`

Let {math}`\mathcal{M}` be a maximal k-subgroup of space group {math}`\mathcal{G}`:
1. {math}`\mathcal{T}(\mathcal{M})` is a maximal {math}`\mathcal{G}`-invariant subgroup of {math}`\mathcal{T}(\mathcal{G})`.
2. {math}`\mathcal{T}(\mathcal{G}) / \mathcal{T}(\mathcal{M})` is an elementary abelian {math}`p`-group for a prime {math}`p`. That is, {math}`\mathcal{T}(\mathcal{G}) / \mathcal{T}(\mathcal{M}) \simeq (\mathbb{Z} / p \mathbb{Z})^{r}` for a non-negative integer {math}`r`.

The above theorem can be proved as follows.
When {math}`\mathcal{M}` is a maximal k-subgroup, {math}`\mathcal{T}(\mathcal{G}) / \mathcal{T}(\mathcal{M})` has no proper subgroup that is invariant under conjugate action of {math}`\mathcal{G}/\mathcal{T}(\mathcal{M})`.
Thus, {math}`\mathcal{T}(\mathcal{G}) / \mathcal{T}(\mathcal{M})` has no charactristic subgroup [^charactristic].
In general, finitely generated abelian group without charastristic subgroups is an elementary abelian {math}`p`-group for a prime {math}`p` [^aut].

[^charactristic]: A charastiristic subgroup of group {math}`G` is a subgroup that is invariant under all authmorphism of {math}`G`.
[^aut]: If {math}`\mathrm{Aut}(G)` acts on {math}`G` transitively, all elements other than identity have the same prime order.

We need to find a translation subgroup {math}`\mathcal{S} < \mathcal{T}(\mathcal{G})` with
- {math}`\mathcal{T}(\mathcal{G}) / \mathcal{S} \simeq \mathbb{Z}_{p^{r}}` for some prime number {math}`p` and integer {math}`r`.
- Let {math}`L(\mathcal{S}) = \mathbf{M}\mathbb{Z}^{n}`. For all {math}`\mathbf{R} \in \mathcal{P}(\mathcal{G})`, {math}`\mathbf{RM}\mathbb{Z}^{n} = \mathbf{M}\mathbb{Z}^{n}`, that is, {math}`\mathbf{M}^{-1}\mathbf{R}\mathbf{M}` is a unimodular matrix.

Let the Smith normal form of {math}`\mathbf{M}` as {math}`\mathbf{M} = \mathrm{diag}(n_{1}, n_{2}, n_{3})`.
The condition that {math}`\mathcal{T}(\mathcal{G}) / \mathcal{S}` has no proper subgroups requires {math}`n_{1} = n_{2} = 1` and {math}`n_{3}` is a power of a prime number, which follow from the Chinese remainder theorem.

## Minimal supercell problem

Reference {cite}`PhysRevB.100.014303` proposed a method to find supercell matrix {math}`\mathbf{H}` that accommodates all given wave vectors {math}`\mathbf{Q}`:
```{math}
    \min_{ \mathbf{H} \in \mathbb{Z}^{d \times d} } \det \mathbf{H}
    \quad \mathrm{s.t.} \quad
    \mathbf{H}^{\top} \mathbf{q} \equiv \mathbf{0}_{d} \, (\mathrm{mod}\, \mathbb{Z}^{d}) \quad (\forall \mathbf{q} \in \mathbf{Q})
```

Ref. {cite}`PhysRevB.92.184301` gave a specific formula for a commensurate supercell for given single wave vector.

## References

```{bibliography}
:filter: docname in docnames
```

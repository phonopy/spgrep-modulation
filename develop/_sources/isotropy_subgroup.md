# Isotropy subgroup

Consider space group {math}`\mathcal{G}` and its representation
```{math}
    g \phi_{j} = \sum_{i} \phi_{i} \Gamma(g)_{ij} \quad (g \in \mathcal{G}).
```
A subgroup {math}`\mathcal{H} (\leq \mathcal{G})` is called **isotropy subgroup** belonging {math}`\Gamma` if the subduced representation {math}`\Gamma \downarrow \mathcal{H}` has an identity representation [^isotropy_subgroup] {cite}`Howard:ta0003,PhysRevB.43.11010,doi:10.1142/0751`.

When we consider a linear combination of basis functions as {math}`\mathbf{\eta} = \sum_{i} \eta_{i} \phi_{i}`, we call {math}`\mathbf{\eta} = \{ \eta_{i} \}_{i}` as **order parameters**.
Two order parameters {math}`\mathbf{\eta}` and {math}`\mathbf{\eta}'` are equivalent if some operation {math}`g \in \mathcal{G}` exists such that
```{math}
    \mathbf{\eta}' = \mathbf{\Gamma}(g) \mathbf{\eta}.
```

[^isotropy_subgroup]: Consider group {math}`G` acting on space {math}`X`.
    In general, isotropy subgroup is defined as stabilizer of point {math}`x` in {math}`X` as
    ```{math}
        G_{x} = \left\{ g \in G \mid g x = x \right\}.
    ```

## Algorithm to generate isotropy subgroups of space group {cite}`Stokes:vk5013`

For a small representation {math}`\Gamma^{\mathbf{k}\alpha}` of space group {math}`\mathcal{G}`, its maximal isotropy subgroups {math}`\mathcal{S}` are enumerated as follows.

### Determine sublattice {math}`\mathcal{L}_{\mathcal{S}}` and translational subgroup {math}`\mathcal{T}(\mathcal{S})`

The requirements that {math}`\mathcal{S}` is an isotropy subgroup of {math}`\mathcal{G}` is
```{math}
    \exists \mathbf{\eta} \neq \mathbf{0} \, s.t. \,
    \mathbf{\Gamma}^{\mathbf{k}\alpha}(g) \mathbf{\eta} = \mathbf{\eta} \quad (\forall g \in \mathcal{S}).
```

In particular, if translation {math}`(\mathbf{E}, \mathbf{t})` belongs to {math}`\mathcal{S}`, {math}`e^{-i\mathbf{k}\cdot\mathbf{t}} = 1`.
Thus, a translational subgroup of {math}`\mathcal{S}` should be
```{math}
    \mathcal{T}(\mathcal{S}) := \{ (\mathbf{E}, \mathbf{t}) \in \mathcal{G} \mid \mathbf{k} \cdot \mathbf{t} \in 2\pi \mathbb{Z} \}.
```
For later use, let {math}`\mathcal{L}_{\mathcal{S}}` be a sublattice {math}`\mathcal{L}_{\mathcal{S}}` formed by translation parts in {math}`\mathcal{T}(\mathcal{S})`.
Note that, although subgroups of {math}`\mathcal{T}(\mathcal{S})` also satisfy the requirements, there is no need to consider such subgroups because these subgroups show at lower-symmetry {math}`\mathbf{k}` vector.

There are two basis vectors for {math}`\mathcal{T}(\mathcal{S})` that are orthogonal to {math}`\mathbf{k}`.
Let {math}`l` be LCM of denominators of {math}`\frac{1}{2\pi}\mathbf{k}`.
Then, we write the elements of {math}`\mathbf{k}` as
```{math}
    \mathbf{k} = 2\pi \left( \frac{g a_{1}}{l} \frac{g a_{2}}{l} \frac{g a_{3}}{l} \right)^{\top},
```
where {math}`\mathrm{GCD}(a_{1}, a_{2}, a_{3}) = 1`, {math}`\mathrm{GCD}(g, l) = 1` and {math}`1 \leq g < l`.

```{math}
    &\mathbf{k} \cdot \mathbf{t} \in 2\pi \mathbb{Z} \\
    &\Leftrightarrow g (a_{1} \,a_{2} \,a_{3}) \mathbf{t} \equiv 0 \quad (\mathrm{mod}\, l) \\
    &\Leftrightarrow (a_{1} \,a_{2} \,a_{3}) \mathbf{t} \equiv 0 \quad (\mathrm{mod}\, l)
        \quad (\because \mathrm{GCD}(g, l) = 1) \\
```

By solving {math}`(a_{1} \,a_{2} \,a_{3}) \mathbf{t} = l`, we obtain one special solution {math}`(a_{1} \,a_{2} \,a_{3}) \mathbf{t}_{0} = l` and two general solutions {math}`(a_{1} \,a_{2} \,a_{3}) \mathbf{t}_{i} = 0 \, (i=1,2)`.
Then, {math}`\{ n\mathbf{t}_{0}, \mathbf{t}_{1}, \mathbf{t}_{2} \}` spans a lattice
```{math}
    \{ \mathbf{t} \in \mathbb{Z}^{3} \mid (a_{1} \,a_{2} \,a_{3}) \mathbf{t} = nl \}
```
for {math}`n \in \mathbb{Z}`.
Thus, {math}`\{ \mathbf{t}_{0}, \mathbf{t}_{1}, \mathbf{t}_{2} \}` is basis of {math}`\mathcal{L}_{S}`.
By stacking these basis vectors, we can find a transformation matrix {math}`\mathbf{M}`,
```{math}
    \mathcal{L}_{\mathcal{S}} = \{ \mathbf{Mt} \mid (\mathbf{E}, \mathbf{t}) \in \mathcal{G} \}.
```

### Enumerate point group {math}`\mathcal{P}(\mathcal{S})`

Next, we find a subgroup of {math}`\mathcal{P}(\mathcal{G})` which preserve the sublattice {math}`\mathcal{L}_{\mathcal{S}}`,
```{math}
    \mathcal{B}(\mathcal{S})
        &:= \{ \mathbf{R} \in \mathcal{P}(\mathcal{G}) \mid \mathbf{R} \mathcal{L}_{\mathcal{S}} =\mathcal{L}_{\mathcal{S}} \} \\
        &= \{ \mathbf{R} \in \mathcal{P}(\mathcal{G}) \mid \mathbf{M}^{-1}\mathbf{R}\mathbf{M} \,\mbox{is unimodular} \} \\
```

Point group of isotropy subgroup {math}`\mathcal{S}` should be a subgroup of {math}`\mathcal{B}(\mathcal{S})`.

### Enumerate isotropy subgroup {math}`\mathcal{S}`

For given point group {math}`\mathcal{P}(\mathcal{S})` and translational subgroup {math}`\mathcal{T}(\mathcal{S})`, consider the following set
```{math}
    \mathcal{S}
        := \{
            ( \mathbf{R}_{i}, \mathbf{\tau}_{i} + \mathbf{c}_{i} + \mathbf{l} )
            \mid
            \mathbf{R}_{i} \in \mathcal{P}(\mathcal{S}), \mathbf{l} \in \mathcal{L}_{\mathcal{S}}
        \}.
```
Here {math}`\mathbf{c}_{i} \in \mathcal{T}(\mathcal{G})` can be freely chosen.

The condition that {math}`\mathcal{S}` is a subgroup of {math}`\mathcal{G}` is as follows:
```{math}
    &\forall ( \mathbf{R}_{i}, \mathbf{\tau}_{i} + \mathbf{Mt} ), ( \mathbf{R}_{j}, \mathbf{\tau}_{j} + \mathbf{Mt}' ) \in \mathcal{S},
        ( \mathbf{R}_{i}, \mathbf{\tau}_{i} + \mathbf{Mt} )^{-1} ( \mathbf{R}_{j}, \mathbf{\tau}_{j} + \mathbf{Mt}' ) \in \mathcal{S} \\
    &\Leftrightarrow
        \forall ( \mathbf{R}_{i}, \mathbf{\tau}_{i} + \mathbf{Mt} ), ( \mathbf{R}_{j}, \mathbf{\tau}_{j} + \mathbf{Mt}' ) \in \mathcal{S},
        \mathbf{\tau}_{j} + \mathbf{Mt}' - \mathbf{R}_{i}^{-1}(\mathbf{\tau}_{i} + \mathbf{Mt}) \in \mathcal{L}_{\mathcal{S}} \\
    &\Leftrightarrow
        \forall \mathbf{R}_{i}, \mathbf{R}_{j} \in \mathcal{S},
        \exists k \,s.t.\, \mathbf{R}_{i}^{-1} \mathbf{R}_{j} = \mathbf{R}_{k},
        \mathbf{\tau}_{j} - \mathbf{R}_{i}^{-1}\mathbf{\tau}_{i} - \mathbf{\tau}_{k} \in \mathcal{L}_{\mathcal{S}} \\
```

### Determine order-parameter direction

Non-zero order-parameter directions correspond to eigenvectors with eigenvalue 1 of the projection operator,
```{math}
    \frac{1}{|\overline{\mathcal{S}}|} \sum_{ g \in \overline{\mathcal{S}} } \Gamma^{\mathbf{k}\alpha}(g).
```

## References

```{bibliography}
:filter: docname in docnames
```

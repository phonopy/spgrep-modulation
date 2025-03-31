# Band connectivity

Refs. {cite}`PhysRevB.97.035138,PhysRevE.96.023310`

We consider to "sew" a band structure between neighboring {math}`\mathbf{k}`-vectors.

## Problem statement

Let {math}`\mathbf{G}` be a space group.
We call the largest continuous submanifold of reciprocal space that the little co-group of every point in the manifold is isomorphic to each other as **{math}`\mathbf{k}`-manifold**.
We say that two {math}`\mathbf{k}`-manifolds are **connected** if the two {math}`\mathbf{k}`-manifolds have a common {math}`\mathbf{k}`-vector with some specific free parameters [^diff].
We call a manifold of {math}`\mathbf{k}`-vector is of **maximal symmetry** if its little co-group is not a proper subgroup of {math}`\mathbf{k}`-vectors in the connected manifolds.
We refer to a {math}`\mathbf{k}`-vector in {math}`\mathbf{k}`-manifold of maximal symmetry as **maximal {math}`\mathbf{k}`-vector**.

[^diff]: For our aim, this definition does not consider reciprocal lattice translations unlike Ref. {cite}`PhysRevB.97.035138`.

For a {math}`\mathbf{k}`-vector, we consider a pair of irrep formed by eigenvectors and corresponding eigenvalue as a node of a graph:
```{math}
    \mathcal{N}_{\mathbf{k}} = \{ (\Gamma^{\mathbf{k}\alpha}, \lambda) \mid \mbox{some eigenvectors with eigenvalue $\lambda$ form irrep $\Gamma^{\mathbf{k}\alpha}$} \}.
```
Given two sets of nodes {math}`\mathcal{N}_{\mathbf{k}_{1}}` and {math}`\mathcal{N}_{\mathbf{k}_{2}}` from two connected {math}`\mathbf{k}`-manifolds, these nodes are linked via {math}`\mathcal{N}_{\mathbf{k}_{p}}` where a {math}`\mathbf{k}`-manifold containing {math}`\mathbf{k}_{p}` is connected to both {math}`\mathbf{k}_{1}` and {math}`\mathbf{k}_{2}`.

Two nodes {math}`(\Gamma^{\mathbf{k}_{i}\alpha}, \lambda) \in \mathcal{N}_{\mathbf{k}_{i}}` and {math}`(\Gamma^{\mathbf{k}_{p}\beta}, \lambda') \in \mathcal{N}_{\mathbf{k}_{p}}` are connected if the subduced representation of {math}`\Gamma^{\mathbf{k}_{i}\alpha}` contains {math}`\Gamma^{\mathbf{k}_{p}\beta}`,
```{math}
    \Gamma^{\mathbf{k}_{i}\alpha} \downarrow \mathcal{G}^{\mathbf{k}_{p}}
        \simeq \bigoplus_{\beta} m_{\beta} \Gamma^{\mathbf{k}_{\mathbf{p}}\beta}
        \quad (i=1,2).
```

## References

```{bibliography}
:filter: docname in docnames
```

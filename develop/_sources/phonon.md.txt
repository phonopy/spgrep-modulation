# Lattice vibration

## Harmonic phonon

Hamiltonian and commutation relations:
```{math}
  H = \sum_{ l \kappa \mu } \frac{ p_{\mu}(l\kappa)^{2} }{ 2 M_{\kappa} }
    + \frac{1}{2} \sum_{ l \kappa \mu } \sum_{ l' \kappa' \mu' } \Phi_{ \mu \mu' }(l\kappa, l'\kappa') u_{\mu}(l\kappa) u_{\mu'}(l'\kappa')
```

Dynamical matrix [^dynamical_matrix]
```{math}
  D_{\mu \mu'}(\kappa\kappa'; \mathbf{q})
    = \frac{1}{ \sqrt{ M_{\kappa} M_{\kappa'} } } \sum_{l'} \Phi_{ \mu \mu' }(0\kappa, l'\kappa') e^{i \mathbf{q} \cdot ( \mathbf{r}(l'\kappa') - \mathbf{r}(0\kappa) )}
```
We denote {math}`[\mathbf{D}(\mathbf{q})]_{\kappa\mu, \kappa'\mu'} = D_{\mu \mu'}(\kappa\kappa'; \mathbf{q})`, then
```{math}
  \mathbf{D}(-\mathbf{q}) &= \mathbf{D}(\mathbf{q})^{\ast} \\
  \mathbf{D}(\mathbf{q})^{\dagger} &= \mathbf{D}(\mathbf{q}) \quad \mbox{(Hermite)}.
```
Here we use {math}`\Phi_{ \mu' \mu }(0\kappa', -l\kappa) = \Phi_{ \mu \mu' }(0\kappa, l'\kappa')`.

Let normalized eigenvector of {math}`\mathbf{D}(\mathbf{q})` as {math}`[\mathbf{e}(\mathbf{q}\nu)]_{\kappa\mu} = e_{\mu}(\kappa; \mathbf{q}\nu)` with
```{math}
  \mathbf{D}(\mathbf{q}) \mathbf{e}(\mathbf{q}\nu)
    &= \omega_{\mathbf{q}\nu}^{2} \mathbf{e}(\mathbf{q}\nu)
    \quad (\nu = 1, \dots, 3N) \\
  \omega_{-\mathbf{q}\nu}^{2} &= \omega_{\mathbf{q}\nu}^{2}.
```
We can choose as
```{math}
  e_{\mu}(\kappa; -\mathbf{q}\nu) = e_{\mu}(\kappa; \mathbf{q}\nu)^{\ast}.
```
Later we denote {math}`q = (\mathbf{q}, \nu)` and {math}`-q = (-\mathbf{q}, \nu)`.

[^dynamical_matrix]: This is the same phase convention with [phonopy](https://phonopy.github.io/phonopy/formulation.html#dynamical-matrix).
    There is the other formulation for defining dynamical matrix as
    ```{math}
    \Psi_{\mu \mu'}(\kappa\kappa'; \mathbf{q})
      = \frac{1}{ \sqrt{ M_{\kappa} M_{\kappa'} } } \sum_{l'} \Phi_{ \mu \mu' }(0\kappa, l'\kappa') e^{i \mathbf{q} \cdot \mathbf{r}(l')}.
    ```

## Normal coordinates

```{math}
  u_{\mu}(l\kappa)
    &=: \frac{1}{\sqrt{ L^{3} M_{\kappa} }} \sum_{q} Q_{q} e_{\mu}(\kappa; q) e^{i \mathbf{q} \cdot \mathbf{r}(l\kappa)} \\
  p_{\mu}(l\kappa)
    &=: \sqrt{\frac{M_{\kappa}}{L^{3}}} \sum_{q} P_{q} e_{\mu}(\kappa; -q) e^{-i \mathbf{q} \cdot \mathbf{r}(l\kappa)} \\
  Q_{q}
    &= \sum_{ l\kappa\mu } \sqrt{ \frac{M_{\kappa}}{L^{3}} } u_{\mu}(l\kappa) e_{\mu}(\kappa; -q) e^{-i \mathbf{q} \cdot \mathbf{r}(l\kappa)} \\
  P_{q}
    &= \sum_{ l\kappa\mu } \frac{1}{\sqrt{ L^{3}M_{\kappa} }} p_{\mu}(l\kappa) e_{\mu}(\kappa; q) e^{i \mathbf{q} \cdot \mathbf{r}(l\kappa)} \\
  Q_{-q} &= Q_{q}^{\ast} \\
  P_{-q} &= P_{q}^{\ast} \\
  H &= \frac{1}{2} \sum_{q} \left( P_{q} P_{-q} + \omega_{q}^{2} Q_{q} Q_{-q} \right) \\
```

Creation and annihilation operators
```{math}
  A_{q}
    &:= \frac{1}{\sqrt{2\hbar\omega_{q}}} \left( \omega_{q} Q_{q} + i P_{-q} \right) \\
  A_{q}^{\dagger}
    &:= \frac{1}{\sqrt{2\hbar\omega_{q}}} \left( \omega_{q} Q_{-q} - i P_{q} \right) \\
  Q_{q}
    &= \sqrt{\frac{\hbar}{2\omega_{q}}} \left( A_{q} + A_{-q}^{\dagger} \right) \\
  P_{q}
    &= i \sqrt{\frac{\hbar\omega_{q}}{2}} \left( A_{-q}^{\dagger} - A_{q} \right) \\
  \left[ A_{q}, A_{q'}^{\dagger} \right] &= \delta_{qq'} \\
  \left[ A_{q}, A_{q'} \right] &= 0 \\
  \left[ A_{q}^{\dagger}, A_{q'}^{\dagger} \right] &= 0 \\
  H &= \frac{1}{2} \sum_{q} \hbar \omega_{q} \left( A_{q}^{\dagger}A_{q} + \frac{1}{2} \right)
```

In canonical ensemble:
```{math}
  \langle A_{q}^{\dagger}A_{q} \rangle_{\beta}
    &= \frac{1}{ e^{\beta \hbar \omega_{q}} - 1 } \\
  \langle |Q_{q}|^{2} \rangle_{\beta}
    &= \frac{\hbar}{2\omega_{q}} \left( 2 \langle A_{q}^{\dagger}A_{q} \rangle_{\beta} + 1 \right)
```


## Action on displacements

We define left group action for positions {math}`\mathbf{r}(l\kappa) = \mathbf{r}(l) + \mathbf{r}(0\kappa)` by {math}`g = (\mathbf{R}_{g}, \mathbf{\tau}_{g}) \in \mathcal{G}` as
```{math}
  g \mathbf{r}(l\kappa)
    := \mathbf{R}_{g} \mathbf{r}(l\kappa) + \mathbf{\tau}_{g}.
```
We denote that site {math}`\mathbf{r}(0, \kappa)` is transformed to {math}`\mathbf{r}(0, g \kappa ) + \mathbf{h}_{g}(\kappa)` by symmetry operation {math}`g`.
Then,
```{math}
    g \mathbf{r}(\mathbf{l}, \kappa) &= \mathbf{R}_{g} \mathbf{r}(l) + \mathbf{r}(0, g\kappa) + \mathbf{h}_{g}(\kappa) \\
```
% \mathbf{h}_{g^{-1}}(\kappa) &= - \mathbf{p}_{g}^{-1} \mathbf{h}_{g}(g^{-1}\kappa)

We define left group action for displacement at {math}`\mathbf{r}(l\kappa)` as [^displacement_action]
```{math}
  g u_{\mu}(\mathbf{r}(l\kappa))
    := \sum_{\nu} [\mathbf{R}_{g}]_{\mu\nu} u_{\nu}(g \mathbf{r}(l\kappa)).
```

[^displacement_action]: This definition actually satisfies the condition of left group action
    ```{math}
        \left[ g \left( g' \mathbf{u}(\mathbf{r}(l\kappa)) \right) \right]_{\mu}
        &= \left[ g \left\{ \sum_{\nu} R_{g', \mu'\nu} u_{\nu}( g'\mathbf{r}(l\kappa) ) \right\}_{\mu'} \right]_{\mu} \\
        &= \sum_{ \mu'\nu } R_{g, \mu\mu'} R_{g', \mu'\nu} u_{\nu}( gg'\mathbf{r}(l\kappa) ) \\
        &= \sum_{ \nu } R_{gg', \mu\nu} u_{\nu}( gg'\mathbf{r}(l\kappa) ) \\
        &= \left[ (gg') \mathbf{u}(\mathbf{r}(l\kappa)) \right]_{\mu}.
    ```

Consider Fourier transformation of {math}`\mathbf{u}(\mathbf{r}(l\kappa))`
```{math}
  \mathbf{u}(\kappa; \mathbf{q})
    &:= \sqrt{\frac{M_{\kappa}}{L^{3}}} \sum_{l} \mathbf{u}(\mathbf{r}(l\kappa)) e^{ i \mathbf{q} \cdot \mathbf{r}(l) } \\
  \mathbf{u}(\mathbf{r}(l\kappa))
    &= \frac{1}{\sqrt{M_{\kappa}L^{3}}} \sum_{\mathbf{q}} \mathbf{u}(\kappa; \mathbf{q}) e^{ -i \mathbf{q} \cdot \mathbf{r}(l) } \\
  g u_{\mu}(\kappa; \mathbf{q})
    &= \sqrt{\frac{M_{\kappa}}{L^{3}}} \sum_{l} g u_{\mu}(\mathbf{r}(l\kappa)) e^{i \mathbf{q}\cdot \mathbf{r}(l)} \\
    &= \sqrt{\frac{M_{\kappa}}{L^{3}}} \sum_{l}\sum_{\nu}
      R_{g,\mu\nu} u_{\nu}\left( \mathbf{R}_{g}\mathbf{r}(l) + \mathbf{r}(0, g\kappa) + \mathbf{h}_{g}(\kappa) \right)
      e^{i \mathbf{q}\cdot \mathbf{r}(l)} \\
      &\quad (\mathbf{r}(l') := \mathbf{R}_{g}\mathbf{r}(l) + \mathbf{h}_{g}(\kappa) ) \\
    &= \sqrt{\frac{M_{\kappa}}{L^{3}}} \sum_{l'}\sum_{\nu} R_{g,\mu\nu} u_{\nu}\left( \mathbf{r}(l', g\kappa) \right) e^{i \mathbf{q}\cdot \mathbf{R}_{g}^{-1}(\mathbf{r}(l') - \mathbf{h}_{\kappa} ) } \\
    &= \sqrt{\frac{M_{\kappa}}{L^{3}}} \sum_{l'}\sum_{\nu} R_{g,\mu\nu} u_{\nu}\left( \mathbf{r}(l', g\kappa) \right) e^{i \mathbf{R}_{g}\mathbf{q}\cdot (\mathbf{r}(l') - \mathbf{h}_{\kappa} ) }
        \quad (\because \mathbf{R}_{g} \in O(3)) \\
    &= \sum_{\kappa'\mu'} u_{\mu'}(\kappa'; \mathbf{R}_{g} \mathbf{q} ) \Gamma_{\kappa'\mu'; \kappa\mu}^{\mathbf{q}}(g)
      \quad (\because M_{g\kappa} = M_{\kappa}),
```
where
```{math}
:label: dynamical_matrix_rep
  \Gamma_{\kappa'\mu'; \kappa\mu}^{\mathbf{q}}(g)
    &:= \exp \left( -i \mathbf{R}_{g} \mathbf{q} \cdot \mathbf{h}_{g}(\kappa) \right) [\mathbf{R}_{g}]_{\mu'\mu} \delta_{ g\kappa, \kappa' } \\
    &:= \exp \left( -2 \pi i \mathbf{R}_{g, f}^{\top} \mathbf{q}_{f} \cdot \mathbf{h}_{g, f}(\kappa) \right) [\mathbf{R}_{g}]_{\mu'\mu} \delta_{ g\kappa, \kappa' } \\
    &\quad (
      \mathbf{R}_{g, f} := \mathbf{A}^{-1} \mathbf{R}_{g} \mathbf{A},
      \mathbf{q} =: 2\pi \mathbf{A}^{-\top} \mathbf{q}_{f},
      \mathbf{v} =: \mathbf{A} \mathbf{v}_{f}
    )
```
Equation {eq}`dynamical_matrix_rep` is essentially the same with Eq. (2.37) of {cite}`RevModPhys.40.1`.

We write matrix representation {math}`[\mathbf{\Gamma}^{\mathbf{q}}(g)]_{ \kappa'\mu'; \kappa\mu } := \Gamma_{\kappa'\mu'; \kappa\mu}^{\mathbf{q}}(g)`.
Then,
```{math}
    \left[ \mathbf{\Gamma}^{ \mathbf{q}}(gg') \right]_{ \kappa'\mu', \kappa\mu }
    &= \delta_{gg'\kappa, \kappa'} R_{gg', \mu'\mu} \exp \left( -i \mathbf{R}_{gg'}\mathbf{q} \cdot \mathbf{h}_{gg'}(\kappa) \right) \\
    &= \sum_{ \kappa''\nu }
        \delta_{g\kappa'', \kappa'} \delta_{g'\kappa, \kappa''}
        R_{g, \mu'\nu} R_{g', \nu\mu}
        \exp \left( -i \mathbf{R}_{g'}\mathbf{q} \cdot \mathbf{h}_{g'}(\kappa) \right)
        \exp \left( -i \mathbf{R}_{g}\mathbf{R}_{g'}\mathbf{q} \cdot \mathbf{h}_{g}(g\kappa') \right) \\
    &= \left[ \mathbf{\Gamma}^{ \mathbf{R}_{g'} \mathbf{q}}(g) \mathbf{\Gamma}^{ \mathbf{q}}(g') \right]_{ \kappa'\mu', \kappa\mu } \\
  \mathbf{\Gamma}^{\mathbf{q}}(g)^{\dagger} \mathbf{\Gamma}^{ \mathbf{q}}(g)
    &= \mathbf{1} \quad \mbox{(Unitary)} \\
  \Gamma^{\mathbf{q}}((E, \mathbf{t}))_{ \kappa'\mu', \kappa\mu }
    &= \exp \left( -i \mathbf{q} \cdot \mathbf{t} \right) \delta_{\mu'\mu} \delta_{ \kappa, \kappa' }.
```

Fourier transformation of force constants
```{math}
  \Phi_{\mu\mu'}(\kappa\kappa'; \mathbf{q})
    &:= \frac{1}{\sqrt{M_{\kappa}M_{\kappa'}}} \sum_{l'} \Phi_{\mu\mu'}(0\kappa; l'\kappa') e^{ i \mathbf{q} \cdot \mathbf{r}(l') } \\
  \sum_{ l l' } \Phi_{ \mu \mu' }(l\kappa, l'\kappa') u_{\mu}(l\kappa) u_{\mu'}(l'\kappa')
    &= \sum_{\mathbf{q}} \Phi_{\mu\mu'}(\kappa\kappa'; \mathbf{q}) u_{\mu}(\kappa; \mathbf{q}) u_{\mu'}(\kappa'; -\mathbf{q}) \\
```

The condition that potential energy is invariant under symmetry operations is rewritten as [^fourier_force_constant]
```{math}
  \mathbf{\Phi}(\mathbf{R}_{g} \mathbf{q})
    = \mathbf{\Gamma}^{\mathbf{q}}(g) \mathbf{\Phi}(\mathbf{q}) \mathbf{\Gamma}^{\mathbf{q}}(g)^{\dagger}.
```

[^fourier_force_constant]: The derivation is as follows:
    ```{math}
    \sum_{ l\kappa\mu l'\kappa'\mu' } \Phi_{ \mu \mu' }(l\kappa, l'\kappa') u_{\mu}(l\kappa) u_{\mu'}(l'\kappa')
    &= \sum_{ l\kappa\mu l'\kappa'\mu' } \Phi_{ \mu \mu' }(l\kappa, l'\kappa') gu_{\mu}(l\kappa) gu_{\mu'}(l'\kappa') \\
    &= \dots = \sum_{  } \left[ \mathbf{\Gamma}^{\mathbf{q}}(g) \mathbf{\Phi}(\mathbf{q}) \mathbf{\Gamma}^{\mathbf{q}}(g)^{\dagger} \right]_{\kappa\mu, \kappa'\mu'} u_{\mu}(\kappa; \mathbf{R}_{g}\mathbf{q}) u_{\mu'}(\kappa'; -\mathbf{R}_{g}\mathbf{q}).
    ```

## Small representation of {math}`\mathcal{G}^{\mathbf{q}}`

For {math}`h, h' \in \mathcal{G}^{\mathbf{q}}`,
```{math}
  \mathbf{\Gamma}^{ \mathbf{q}}(h) \mathbf{\Gamma}^{ \mathbf{q}}(h')
    &= \mathbf{\Gamma}^{ \mathbf{q}}(hh') \\
  \mathbf{\Phi}(\mathbf{q})
    &= \mathbf{\Gamma}^{\mathbf{q}}(h) \mathbf{\Phi}(\mathbf{q}) \mathbf{\Gamma}^{\mathbf{q}}(h)^{\dagger}
    \quad (\forall h \in \mathcal{G}^{\mathbf{q}}).
```

We can introduce projective representation {math}`\overline{\Gamma}^{ \mathbf{q}}`,
```{math}
  \mathbf{\Gamma}^{ \mathbf{q}}(h)
    &=: e^{ -i \mathbf{q} \cdot \mathbf{v}_{h} } \overline{\mathbf{\Gamma}}^{ \mathbf{q}}(h) \\
  \overline{\mathbf{\Gamma}}^{ \mathbf{q}}(h) \overline{\mathbf{\Gamma}}^{ \mathbf{q}}(h')
    &= e^{ -i \mathbf{q} \cdot ( \mathbf{R}_{h} \mathbf{v}_{h'} - \mathbf{v}_{h'} ) } \overline{\mathbf{\Gamma}}^{ \mathbf{q}}(hh') \\
  \overline{\mathbf{\Gamma}}^{ \mathbf{q}}((E, \mathbf{t}))
    &= \mathbf{1}
```
Thus, we only need to consider projective representation {math}`\overline{\mathbf{\Gamma}}^{ \mathbf{q}}` for little co-group {math}`\overline{\mathcal{G}}^{\mathbf{q}} \simeq \mathcal{G}^{\mathbf{q}} / \mathcal{T}`.
The decomposition of the projective representation
```{math}
  \overline{\Gamma}^{\mathbf{q}} &= \sum_{\alpha} \sum_{\sigma} \overline{\Gamma}^{\mathbf{q}\alpha\sigma} \\
```
can be performed with [spgrep](https://github.com/spglib/spgrep), where {math}`\alpha` represent irrep and {math}`\sigma = 1,\dots, m_{\alpha}` distinguish equivalent irreps to {math}`\alpha`.
The corresponding small representation of {math}`\mathcal{G}^{\mathbf{q}}` is obtained by {math}`\mathbf{\Gamma}^{ \mathbf{q}\omega}(h) := e^{ -i \mathbf{q} \cdot \mathbf{v}_{h} } \overline{\mathbf{\Gamma}}^{ \mathbf{q}\omega }(h)`.
We call orthonormal basis vectors {math}`f_{\mu}(\kappa; \mathbf{q}\alpha\sigma\nu)` forming irrep {math}`\Gamma^{\mathbf{q}\alpha}` as *modified eigenvectors*:
```{math}
  h f_{\mu}(\kappa; \mathbf{q}\alpha\sigma\nu)
    &= \sum_{\nu'} f_{\mu}(\kappa; \mathbf{q}\alpha\sigma\nu') \Gamma^{\mathbf{q}\alpha}(h)_{\nu',\nu}
    \quad (h \in \mathcal{G}^{\mathbf{q}}, \nu = 1, \dots, d_{\alpha}) \\
  \sum_{\kappa\mu} f_{\mu}(\kappa; \mathbf{q}\alpha\sigma\nu)^{\ast} f_{\mu}(\kappa; \mathbf{q}\alpha\sigma\nu')
    &= \delta_{\nu\nu'}
    \quad (\nu, \nu' = 1, \dots, d_{\alpha})
```

We can subdivide eigenvectors further by decomposing {math}`\Gamma^{\mathbf{q}}` into irreps,
```{math}
  \left[ F^{\mathbf{q}\alpha\sigma} \right]_{\kappa\mu, \nu}
    &:= f_{\mu}(\kappa; \mathbf{q}\alpha\sigma\nu) \\
  \mathbf{F}^{\mathbf{q}\alpha\sigma \dagger} \mathbf{\Gamma}^{\mathbf{q}}(h) \mathbf{F}^{\mathbf{q}\alpha\sigma}
    &= \mathbf{\Gamma}^{\mathbf{q}\alpha}(h)
    \quad (h \in \mathcal{G}^{\mathbf{q}}) \\
```

## Solve dynamical matrix w.r.t. modified eigenvectors

Block-diagonalize fourier transformed force constants:
```{math}
  \mathbf{F}^{\mathbf{q}\alpha}
    &:= \left( \mathbf{F}^{\mathbf{q}\alpha 1} \dots \mathbf{F}^{\mathbf{q}\alpha m_{\alpha}} \right)
    \quad \in \mathbb{C}^{3N \times m_{\alpha}d_{\alpha}} \\
  \mathbf{\Phi}(\mathbf{q}\alpha)
    &:= \mathbf{F}^{\mathbf{q}\alpha \dagger} \mathbf{\Phi}(\mathbf{q}) \mathbf{F}^{\mathbf{q}\alpha}
    \quad \in \mathbb{C}^{ m_{\alpha}d_{\alpha} \times m_{\alpha}d_{\alpha} }, \\
```
where {math}`\mathbf{\Phi}(\mathbf{q}\alpha)` is hermitian.

Diagonalize {math}`\mathbf{\Phi}(\mathbf{q}\alpha)`
```{math}
  \sum_{\sigma' \nu'} \Phi(\mathbf{q}\alpha)_{\sigma\nu, \sigma'\nu'} c(\mathbf{q} \alpha s\lambda)_{\sigma' \nu'}
    &= \omega_{\alpha s}^{2} c(\mathbf{q} \alpha s\lambda)_{\sigma \nu}
    \quad (s = 1, \dots, m_{\alpha}, \lambda = 1, \dots, d_{\alpha}) \\
  \sum_{\sigma\nu} c(\mathbf{q}\alpha s\lambda)_{\sigma\nu}^{\ast} c(\mathbf{q}\alpha s'\lambda')_{\sigma\nu}
    &= \delta_{s s'} \delta_{\lambda \lambda'} \\
  \sum_{s\lambda} c(\mathbf{q}\alpha s\lambda)_{\sigma\nu}^{\ast} c(\mathbf{q}\alpha s\lambda)_{\sigma'\nu'}
    &= \delta_{\sigma \sigma'} \delta_{\nu \nu'}
    \quad (\sigma = 1, \dots, m_{\alpha}, \nu = 1, \dots, d_{\alpha}) \\
```
where {math}`s` labels real eigenvalues {math}`\omega_{\alpha s}^{2}` and {math}`\lambda` labels degenerated eigenvectors.
Here, **we assume each equivalent irrep {math}`\Gamma^{\mathbf{q}\alpha \sigma}` gives a different eigenvalue**.
Also, we choice eigenvectors {math}`\mathbf{c}(\mathbf{q} \alpha s\lambda) := \{ c(\mathbf{q}\alpha s\lambda)_{\sigma\nu} \}_{\sigma\nu}` are mutually orthogonal even within degenerated eigenvalues.

When irrep {math}`\alpha` appears more than once ({math}`m_{\alpha} > 1`), irrep formed by eigenvectors is no longer same as {math}`\Gamma^{\mathbf{q}\alpha}` {cite}`RevModPhys.40.1`:
```{math}
  \tilde{\mathbf{f}}(\mathbf{q}\alpha s \lambda)
    &:= \sum_{\sigma\nu} \mathbf{f}(\mathbf{q}\alpha \sigma \nu) c(\mathbf{q}\alpha s\lambda)_{\sigma \nu} \\
  h \tilde{\mathbf{f}}(\mathbf{q}\alpha s \lambda)
    &= \sum_{\lambda'} \tilde{\mathbf{f}}(\mathbf{q}\alpha s\lambda') \tilde{\Gamma}^{\mathbf{q}\alpha s}(h)_{\lambda'\lambda},
```
where [^eigmode_rep]
```{math}
  \tilde{\Gamma}^{\mathbf{q}\alpha s}(h)_{\lambda'\lambda}
    := \sum_{\sigma\nu\nu'} c(\mathbf{q}\alpha s\lambda')_{\sigma\nu'}^{\ast} \Gamma^{\mathbf{q}\alpha}(h)_{\nu'\nu} c(\mathbf{q}\alpha s\lambda)_{\sigma\nu}
    \quad (h \in \mathcal{G}^{\mathbf{q}}).
```

[^eigmode_rep]: The derivation is cumbersome (to me).
    Let {math}`\mathbf{C}^{\mathbf{q}\alpha s} = ( \mathbf{c}(\mathbf{q}\alpha s 1), \dots, \mathbf{c}(\mathbf{q}\alpha s d_{\alpha}) )`.
    ```{math}
    \tilde{\mathbf{\Gamma}}^{\alpha s}(h)
      &:= \mathbf{C}^{\mathbf{q}\alpha s}
         \begin{pmatrix}
          \mathbf{\Gamma}^{\mathbf{q}\alpha}(h) & & \\
          & \ddots & \\
          & & \mathbf{\Gamma}^{\mathbf{q}\alpha}(h) \\
         \end{pmatrix}
         \mathbf{C}^{\mathbf{q}\alpha s \dagger} \\
    \tilde{\Gamma}^{\alpha s}(h)_{\lambda' \lambda}
      &= \sum_{\sigma' \nu' \sigma \nu}
            c(\mathbf{q}\alpha s \lambda')_{\sigma'\nu'}^{\ast} \delta_{\sigma \sigma'}
            \Gamma^{\mathbf{q}\alpha}(h)_{\nu'\nu}
            c(\mathbf{q}\alpha s \lambda)_{\sigma\nu} \\
      &= \sum_{\sigma \nu' \nu}
            c(\mathbf{q}\alpha s \lambda')_{\sigma\nu'}^{\ast}
            \Gamma^{\mathbf{q}\alpha}(h)_{\nu'\nu}
            c(\mathbf{q}\alpha s \lambda)_{\sigma\nu} \\
    ```

Now go back to the other convention of dynamical matrix:
```{math}
  e_{\mu}(\kappa; \mathbf{q}\alpha s \lambda)
    &:= e^{ -i\mathbf{q} \cdot \mathbf{r}(0\kappa) } \tilde{f}_{\mu}(\kappa; \mathbf{q}\alpha s \lambda) \\
  D_{\mu\mu'}(\kappa\kappa'; \mathbf{q})
    &= e^{ i \mathbf{q} \cdot \left( \mathbf{r}(0\kappa') - \mathbf{r}(0\kappa) \right) } \Phi_{\mu\mu'}(\kappa\kappa'; \mathbf{q}) \\
  \sum_{ \kappa'\mu' } D_{\mu\mu'}(\kappa\kappa'; \mathbf{q}) e_{\mu'}(\kappa'; \mathbf{q}\alpha s\lambda)
    &= \omega_{\alpha s}^{2} e_{\mu}(\kappa; \mathbf{q}\alpha s \lambda) \\
  h \mathbf{e}(\mathbf{q}\alpha s \lambda)
    &= \sum_{\lambda'} \mathbf{e}(\mathbf{q}\alpha s'\lambda') \tilde{\Gamma}^{\mathbf{q}\alpha s}(h)_{\lambda'\lambda}
    \quad (h \in \mathcal{G}^{\mathbf{q}}).
```

## Modulation

Modulation associated with qpoint {math}`\mathbf{q}` and frequency {math}`\mathbf{\omega}_{\mathbf{q}}`:

{math}`\mathbf{q} \neq \mathbf{0}` case:
```{math}
  u^{( \mathbf{q} \alpha s)}_{\mu}(l\kappa)
    &= \frac{1}{\sqrt{L^{3}M_{\kappa}}} \sum_{ \lambda }
      \left(
        Q^{ (\mathbf{q} \alpha s) }_{\lambda} e_{\mu}(\kappa; \mathbf{q} \alpha s \lambda) e^{ i \mathbf{q} \cdot \mathbf{r}(l\kappa) }
        + \mathrm{c.c.}
      \right) \\
  Q^{ (\mathbf{q} \alpha s) }_{\lambda} &\in \mathbb{C} \\
```

{math}`\mathbf{q} = \mathbf{0}` case:
```{math}
  u^{( \mathbf{0} \alpha s)}_{\mu}(l\kappa)
    &= \frac{1}{\sqrt{L^{3}M_{\kappa}}} \sum_{ \lambda }
        Q^{ (\mathbf{0} \alpha s) }_{\lambda} e_{\mu}(\kappa; \mathbf{0} \alpha s \lambda ) \\
  Q^{ (\mathbf{0} \alpha s) }_{\lambda} &\in \mathbb{R} \\
```

% When we sample modulations in grid, we fix the norm of {math}`\mathbf{Q}^{ (\mathbf{q} \alpha s) }`, which is equivalent to fix the norm of the following mass-weighted displacements,
% ```{math}
%   \sum_{l \kappa \mu} M_{\kappa} \left| u^{( \mathbf{q} \alpha s)}_{\mu}(l\kappa) \right|^{2}
%     \propto \sum_{\lambda} \left| Q^{ (\mathbf{q} \alpha s) }_{\lambda} \right|^{2}.
% ```

From unitary arbitrariness of {math}`\tilde{\Gamma}^{\mathbf{q}\alpha s}`, we can choose as {math}`\mathrm{Im}\, Q^{ (\mathbf{q} \alpha s) }_{\lambda=1} = 0`.
Thus, sampling of {math}`\mathbf{Q}^{ (\mathbf{q} \alpha s) }` is attributed to sampling points from unit sphere {math}`S^{2 d_{\alpha} - 2}`.

% chain-adapted mode {cite}`Aroyo:js0048`

## References

```{bibliography}
:filter: docname in docnames
```

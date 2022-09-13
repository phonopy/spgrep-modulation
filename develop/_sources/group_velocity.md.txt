# Symmetry condition for group velocity

Ref. {cite}`birman1974quick`

Group action on modified eigenvectors
```{math}
    g f_{\mu}(\kappa; \mathbf{q}\omega_{\mathbf{q}}\nu)
        &= \sum_{\mu'\nu'} f_{\mu'}(\kappa'; \mathbf{q}\omega_{\mathbf{q}}\nu') \Gamma^{ (\mathbf{q} \omega_{\mathbf{q}}) }(g)_{\kappa'\mu'\nu', \kappa\mu\nu}
        \quad (g \in \mathcal{G}^{\mathbf{q}}) \\
    \omega_{\mathbf{q}}^{2}
        &= \sum_{ \kappa \mu }\sum_{ \kappa' \mu' } f_{\mu}(\kappa; \mathbf{q}\omega_{\mathbf{q}}\nu)^{\ast} \Phi_{\mu\mu'}(\kappa\kappa'; \mathbf{q}) f_{\mu'}(\kappa'; \mathbf{q}\omega_{\mathbf{q}}\nu)
```

Group velocity
```{math}
    \nabla_{q_{\alpha}} \omega_{\mathbf{q}}^{2}
        = \sum_{ \kappa \mu }\sum_{ \kappa' \mu' }
            f_{\mu}(\kappa; \mathbf{q}\omega_{\mathbf{q}}\nu)^{\ast}
            \left( \nabla_{q_{\alpha}} \Phi_{\mu\mu'}(\kappa\kappa'; \mathbf{q}) \right)
            f_{\mu'}(\kappa'; \mathbf{q}\omega_{\mathbf{q}}\nu)
```

If the group velocity has scalar component, that is
```{math}
:label: indicator
\left\langle
    \Gamma^{ (\mathbf{q} \omega_{\mathbf{q}}) \ast} \otimes \Gamma^{(\mathbf{\nabla})} \otimes \Gamma^{ (\mathbf{q} \omega_{\mathbf{q}}) }
    \mid
    \Gamma^{ (\mathrm{id}) }
\right\rangle
\geq 1,
```
it should be zero at {math}`\mathbf{q}`.

The gradient operator {math}`\mathbf{\nabla}_{\mathbf{q}}` behaves as covariant vector:
```{math}
    \mathbf{q}' &= \mathbf{R}^{\top}\mathbf{q} \quad (\mathbf{R} \in O(3)) \\
    \nabla_{q'_{\alpha}}
        &= \sum_{\beta} \frac{ \partial q_{\beta} }{ \partial q'_{\alpha} } \nabla_{q_{\beta}}
        = \sum_{\beta} R_{\alpha\beta} \nabla_{q_{\beta}}
```

The left hand side of Eq. {eq}`indicator`:
```{math}
    \sum_{g \in \mathcal{G}^{\mathbf{q}} / \mathcal{T} }
        | \chi^{(\mathbf{q} \omega_{\mathbf{q}})}(g) |^{2} \chi^{(\mathbf{\nabla})}(\mathbf{R}_{g})
```

## References

```{bibliography}
:filter: docname in docnames
```

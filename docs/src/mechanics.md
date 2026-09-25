# Applied knot mechanics

Engineering units throughout: lengths in **mm**, forces in **N**, modulus
and stress in **MPa** (1 GPa = 1000 MPa = 1000 N/mm²).

## The superposition formula

At the strand surface,

```math
\sigma_{\text{total}}(s) = \frac{F(s)}{A} + \frac{E\,\kappa(s)\,d}{2},
\qquad A = \pi (d/2)^2,
```

where ``\kappa = 1/R`` is the centreline curvature. Inside the knot core
the carried tension decays by capstan friction,
``F(s) = F_0 e^{-\mu\theta(s)}`` with ``\theta`` the accumulated contact
angle. This is the theory–practice bridge: the same curvature field the
differential geometry computes drives the failure prediction.

## Rope efficiency

`breaking_force(rod)` solves for the load at which the peak superposed
stress reaches UTS (bending stress is load-independent, axial stress scales
linearly with load); `knot_efficiency` divides by the straight-strand
strength ``\sigma_{\mathrm{uts}} A``. Empirical rope knots sit near
0.60–0.75; the model reproduces that range when the geometry carries a
tight nipping turn.

## Topological integrals

`writhe_of` and `linking_of` evaluate the Gauss double integral by midpoint
quadrature; `calugareanu_check` builds a twisted ribbon and verifies
``\mathrm{Lk} = \mathrm{Tw} + \mathrm{Wr}`` numerically — the identity that
underlies DNA topology and the protistology module.

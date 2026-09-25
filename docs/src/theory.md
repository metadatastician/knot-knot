# Theory layer

## Diagram conventions

KnotKnot uses the **KnotAtlas PD convention** (the same as KnotTheory.jl):
a crossing is `(a, b, c, d, sign)` with arcs listed counter-clockwise; the
under-strand enters at `a` and exits at `c`; the over-strand enters at `d`
and exits at `b`; `sign` is +1 for a positive (right-hand) crossing.

## Invariants

- **Jones** via the Kauffman bracket state sum on PD slots: A-smoothings
  pair slots (1,2)+(3,4), B-smoothings pair (2,3)+(4,1);
  ``V(t) = (-A)^{-3w}\langle D\rangle`` under ``A = t^{-1/4}``.
- **Alexander** via Fox calculus on the Wirtinger presentation; the
  (n−1)×(n−1) minor determinant is computed by fraction-free Bareiss
  elimination over ``\mathbb{Z}[t]`` and canonicalised to a symmetric window
  with ``\Delta(1) = +1`` (Conway normalisation).
- **Conway** from Alexander by ``\Delta(t) = \nabla(t^{1/2} - t^{-1/2})``,
  with the antisymmetric branch for even-component links (odd-span
  ``\Delta``), so the Hopf link gives ``\nabla = z``.
- **Signature** from the symmetrised band-based Seifert matrix;
  **genus** from ``g = (c - s + 1)/2``.

## Quandles

A finite quandle is a Cayley table checked against the three axioms
(idempotence, right-invertibility, self-distributivity). The fundamental
quandle presentation of a diagram has one generator per Wirtinger arc; at a
positive crossing with over-arc `y`, entering under-arc `x` and exiting
under-arc `z`: ``z = x \triangleright y`` (inverse operation at negative
crossings). `coloring_count` counts homomorphisms into a finite quandle —
the ``R_3`` count is classical tricolourability.

## Skein relations

`switch_crossing` and `smooth_crossing` build the L₊/L₋/L₀ triple at any
crossing; `verify_conway_skein` evaluates
``\nabla(L_+) - \nabla(L_-) - z\,\nabla(L_0)`` and the test suite asserts
this vanishes on every crossing of the trefoil and figure-eight diagrams.

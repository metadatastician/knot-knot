<!-- coderabbit-ai-generated:v1 -->

<!-- coderabbit-ai-generated-section:v1:df5ecf71af91f515cf1883a1f1a698a5b2e4d2af41d27693e54fa01f2a068cc3 -->

<!--
// SPDX-License-Identifier: CC-BY-SA-4.0
// Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath)
-->
# Protistological topology

Protists are where knot theory meets cell biology: eukaryotic flagella beat
in helical and planar waveforms whose curvature fields are exactly the
objects the mechanics layer computes, and kinetoplastid protists
(trypanosomes and relatives) pack their mitochondrial genome — the
kinetoplast — into a chain-mail of thousands of catenated DNA minicircles.

## Flagellar waveforms

```julia
using KnotKnot
frame = flagellum(:helical; length = 50.0, radius = 2.0, wavelength = 20.0)
beat_summary(frame)
```

Kinds: `:planar` (sinusoidal, sea-urchin-sperm style), `:helical`
(trypanosome style), `:tipped` (amplitude envelope, a breaststroke proxy).
`beat_summary` reports arclength, curvature statistics, total curvature and
the bending load ``\int \kappa^2 ds`` proportional to the elastic power
a motor ensemble must deliver.

## Kinetoplast catenanes

```julia
using KnotKnot
chain = kinetoplast_chain(4)
linking_matrix(chain)
```

builds alternating-plane minicircle rings in the chainmail spacing window
``r < \text{spacing} < 2r`` — adjacent rings link once, next-nearest
rings link zero — and returns the pairwise Gauss-integral linking matrix,
the topological observable kinetoplast disassembly must resolve.

The module is deliberately self-contained (geometry only) so it can split
out as the future **Proctist.jl** without dragging the engineering layer.


<!-- /coderabbit-ai-generated:v1 -->
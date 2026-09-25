# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Skein-relation machinery on planar diagrams.
#
# The Conway potential is characterised by
#     nabla(L+) - nabla(L-) = z * nabla(L0)
# with nabla(unknot) = 1 and nabla(split link) = 0. `switch_crossing`
# produces L+ / L- by a genuine crossing change (the shadow's cyclic stub
# order is preserved while the over/under pairs swap); `smooth_crossing`
# produces L0 by the oriented resolution, which pairs each entering arc
# with the opposite-strand exit fixed by the diagram's slot orientation.
# Resolved arc labels are merged with union-find.

"""
    switch_crossing(d::PlanarDiagram, index::Int, sign::Int) -> PlanarDiagram

The diagram obtained by forcing crossing `index` to the given sign
(+1 = positive, -1 = negative). This is the L+ / L- side of a skein
triple. When the sign actually changes, the crossing is switched
geometrically: the four stubs keep their cyclic order, the over/under
strand pairs swap, and the new slot-1 arc is whichever of the old
over-arcs enters the crossing (per `slot_orientation`).
"""
function switch_crossing(d::PlanarDiagram, index::Int, sign::Int)::PlanarDiagram
    1 <= index <= length(d.crossings) || error("crossing index out of range")
    sign in (-1, 1) || error("sign must be +1 or -1")
    crossings = Crossing[c for c in d.crossings]
    target = crossings[index]
    if target.sign != sign
        orient = slot_orientation(d)
        base = 4 * (index - 1)
        a, b, c, dd = target.arcs
        new_arcs = orient[base + 2] == 1 ? (b, c, dd, a) : (dd, a, b, c)
        crossings[index] = Crossing(new_arcs, sign)
    end
    PlanarDiagram(crossings, d.components)
end

"""
    smooth_crossing(d::PlanarDiagram, index::Int) -> PlanarDiagram

The oriented smoothing (L0) of crossing `index`: the crossing is deleted
and its four arcs are reconnected by the oriented resolution — the
under-strand entry is paired with the over-strand exit and vice versa, as
fixed by `slot_orientation`. Remaining arc labels are merged with
union-find so strands continue smoothly through the resolved site.
"""
function smooth_crossing(d::PlanarDiagram, index::Int)::PlanarDiagram
    1 <= index <= length(d.crossings) || error("crossing index out of range")
    target = d.crossings[index]
    a, b, c, dd = target.arcs
    orient = slot_orientation(d)
    base = 4 * (index - 1)
    u_in, u_out = orient[base + 1] == 1 ? (a, c) : (c, a)
    o_in, o_out = orient[base + 2] == 1 ? (b, dd) : (dd, b)
    pairs = [(u_in, o_out), (o_in, u_out)]

    all_arcs = arcs_of(d)
    idx = Dict{Int,Int}(x => i for (i, x) in enumerate(all_arcs))
    parent = collect(1:length(all_arcs))
    rnk = zeros(Int, length(all_arcs))
    for (u, v) in pairs
        _uf_union!(parent, rnk, idx[u], idx[v])
    end
    remap = Dict{Int,Int}()
    for x in all_arcs
        remap[x] = all_arcs[_uf_find!(parent, idx[x])]
    end

    crossings = Crossing[]
    for (i, cr) in enumerate(d.crossings)
        i == index && continue
        push!(crossings, Crossing(ntuple(j -> remap[cr.arcs[j]], 4), cr.sign))
    end
    PlanarDiagram(crossings, d.components)
end

"""
    skein_triple(d::PlanarDiagram, index::Int)

The skein triple `(L_plus, L_minus, L_zero)` at crossing `index`.
"""
function skein_triple(d::PlanarDiagram, index::Int)
    (
        switch_crossing(d, index, 1),
        switch_crossing(d, index, -1),
        smooth_crossing(d, index),
    )
end

"""
    verify_conway_skein(d::PlanarDiagram, index::Int) -> LaurentPoly

Evaluate ``\\nabla(L_+) - \\nabla(L_-) - z\\,\\nabla(L_0)`` at crossing
`index`. The Conway relation holds exactly when the result is the zero
polynomial.
"""
function verify_conway_skein(d::PlanarDiagram, index::Int)::LaurentPoly
    lp, lm, l0 = skein_triple(d, index)
    np = conway_polynomial(lp)
    nm = conway_polynomial(lm)
    n0 = conway_polynomial(l0)
    z_n0 = lp_shift(n0, 1)
    lp_trim((np - nm) - z_n0)
end

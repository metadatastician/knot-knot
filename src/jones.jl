# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Jones polynomial via the Kauffman bracket state sum on PD codes.
#
# Bracket rules (crossing slots 1..4 = arcs a,b,c,d counter-clockwise):
#   A-smoothing pairs slots (1,2) and (3,4);
#   B-smoothing pairs slots (2,3) and (4,1).
# Each state contributes A^(#A - #B) * d^(loops - 1) with d = -A^2 - A^-2.
# With this smoothing labelling a positive Reidemeister-I kink multiplies
# the bracket by -A^-3, so the writhe correction is the +w power:
# V(t) = (-A^3)^w <D> evaluated at A = t^(-1/4). For a knot diagram every
# resulting t-exponent is integral. Pinned: the all-positive trefoil PD
# gives -t^-4 + t^-3 + t^-1 (KnotAtlas / KnotTheory.jl convention).

const MAX_BRACKET_CROSSINGS = 20

"""
    bracket_polynomial(d::PlanarDiagram) -> LaurentPoly

The Kauffman bracket ``\\langle D \\rangle`` in the variable ``A``.
Invariant under Reidemeister II and III; a Reidemeister I kink multiplies
it by ``-A^{\\pm 3}``. Exponential in the crossing number (bounded by
`MAX_BRACKET_CROSSINGS`).
"""
function bracket_polynomial(d::PlanarDiagram)::LaurentPoly
    n = length(d.crossings)
    n == 0 && return LaurentPoly(0 => 1)
    n > MAX_BRACKET_CROSSINGS && throw(ArgumentError(
        "bracket state sum limited to $MAX_BRACKET_CROSSINGS crossings (got $n)",
    ))

    pos = arc_positions(d)
    pairs = Tuple{Int,Int}[]
    for p in values(pos)
        length(p) == 2 && push!(pairs, (p[1], p[2]))
    end

    function count_loops(state_pairs::Vector{Tuple{Int,Int}})::Int
        adj = Dict{Int,Vector{Int}}()
        for (x, y) in state_pairs
            push!(get!(adj, x, Int[]), y)
            push!(get!(adj, y, Int[]), x)
        end
        seen = Set{Int}()
        loops = 0
        for node in keys(adj)
            node in seen && continue
            stack = [node]
            while !isempty(stack)
                cur = pop!(stack)
                cur in seen && continue
                push!(seen, cur)
                for nb in get(adj, cur, Int[])
                    nb in seen || push!(stack, nb)
                end
            end
            loops += 1
        end
        loops
    end

    function expand(idx::Int, state_pairs::Vector{Tuple{Int,Int}})::LaurentPoly
        if idx > n
            loops = count_loops(state_pairs)
            poly = LaurentPoly(0 => 1)
            for _ in 1:(loops - 1)
                nxt = LaurentPoly()
                for (e, c) in poly
                    nxt[e + 2] = get(nxt, e + 2, 0) - c
                    nxt[e - 2] = get(nxt, e - 2, 0) - c
                end
                poly = nxt
            end
            return poly
        end
        s0 = 4 * (idx - 1)
        a_pairs = vcat(state_pairs, [(s0 + 1, s0 + 2), (s0 + 3, s0 + 4)])
        b_pairs = vcat(state_pairs, [(s0 + 2, s0 + 3), (s0 + 4, s0 + 1)])
        pa = expand(idx + 1, a_pairs)
        pb = expand(idx + 1, b_pairs)
        r = LaurentPoly()
        for (e, c) in pa
            r[e + 1] = get(r, e + 1, 0) + c
        end
        for (e, c) in pb
            r[e - 1] = get(r, e - 1, 0) + c
        end
        r
    end

    expand(1, pairs)
end

"""
    jones_polynomial(d::PlanarDiagram) -> LaurentPoly

The Jones polynomial ``V(t)`` of a knot diagram, normalised so that
``V(\\mathrm{unknot}) = 1``. Exponents are integral in ``t`` for knots.

Pinned value: the right-hand trefoil gives ``-t^{-4} + t^{-3} + t^{-1}``.
"""
function jones_polynomial(d::PlanarDiagram)::LaurentPoly
    n = length(d.crossings)
    n == 0 && return LaurentPoly(0 => 1)
    br = bracket_polynomial(d)
    w = writhe(d)
    sgn = isodd(w) ? -1 : 1
    jones = LaurentPoly()
    for (e, c) in br
        # V(t) = (-A^3)^w <D> with A = t^(-1/4): A^(e + 3w) -> t^(-(e + 3w)/4).
        texp = -(e + 3 * w)
        texp % 4 == 0 || error("non-integral Jones exponent: input is not a knot diagram")
        q = div(texp, 4)
        jones[q] = get(jones, q, 0) + sgn * c
    end
    lp_trim(jones)
end

# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Braid words and their closures.
#
# Generators: (i, +1) is sigma_i (strand at position i passes OVER the
# strand at position i+1); (i, -1) is its inverse. The closure connects
# bottom strand k back to top strand k, giving a PD code in the same
# KnotAtlas convention used throughout KnotKnot (and by KnotTheory.jl's
# from_braid_word).

"""
    Braid(generators::Vector{Tuple{Int,Int}})

A braid word as `(strand_index, sign)` pairs; `sign` +1 for ``\\sigma_i``
and -1 for ``\\sigma_i^{-1}``.
"""
struct Braid
    generators::Vector{Tuple{Int,Int}}
end

"""
    parse_braid_word(word::AbstractString) -> Braid

Parse TANGLE-style braid words: `s1.s1.s1` (positive generators) and
`S2` for inverses. Whitespace around tokens is tolerated.
"""
function parse_braid_word(word::AbstractString)::Braid
    gens = Tuple{Int,Int}[]
    for raw in split(strip(String(word)), ".")
        token = strip(String(raw))
        isempty(token) && continue
        first_char = token[1:1]
        first_char in ("s", "S") || error("invalid braid generator: '$token'")
        idx = parse(Int, token[2:end])
        idx >= 1 || error("braid generator index must be >= 1: '$token'")
        push!(gens, (idx, first_char == "S" ? -1 : 1))
    end
    Braid(gens)
end

"""
    braid_closure(b::Braid) -> PlanarDiagram

Close the braid (bottom strand k to top strand k) and return the planar
diagram of the closure. The empty braid closes to the unknot.
"""
function braid_closure(b::Braid)::PlanarDiagram
    isempty(b.generators) && return PlanarDiagram(Crossing[], Vector{Vector{Int}}())
    n_strands = maximum(g[1] for g in b.generators) + 1
    current_arc = collect(1:n_strands)
    next_arc = n_strands + 1
    crossings = Crossing[]

    for (i, sgn) in b.generators
        in_i = current_arc[i]
        in_i1 = current_arc[i + 1]
        out_i = next_arc
        out_i1 = next_arc + 1
        next_arc += 2
        if sgn > 0
            # Strand entering at position i goes OVER, exiting at position
            # i+1; strand i+1 passes under (in_i1 -> out_i). The cyclic slot
            # order is the one consistent with the KnotAtlas PD pins used in
            # knot_table.jl (verified: closures of s1^3 and s1^5 reproduce
            # the pinned Jones/Alexander values of 3_1 and 5_1).
            push!(crossings, Crossing((in_i1, in_i, out_i, out_i1), 1))
        else
            # sigma_i^-1: strand i passes under (in_i -> out_i1), strand
            # i+1 goes over (in_i1 -> out_i).
            push!(crossings, Crossing((in_i, out_i, out_i1, in_i1), -1))
        end
        current_arc[i] = out_i
        current_arc[i + 1] = out_i1
    end

    rename = Dict{Int,Int}()
    for k in 1:n_strands
        current_arc[k] != k && (rename[current_arc[k]] = k)
    end
    if !isempty(rename)
        crossings = Crossing[
            Crossing(ntuple(j -> get(rename, c.arcs[j], c.arcs[j]), 4), c.sign) for c in crossings
        ]
    end
    PlanarDiagram(crossings, Vector{Vector{Int}}())
end

function Base.show(io::IO, b::Braid)
    tokens = [sgn > 0 ? "s$i" : "S$i" for (i, sgn) in b.generators]
    print(io, "Braid(", join(tokens, "."), ")")
end

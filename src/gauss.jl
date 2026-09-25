# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Oriented Gauss words derived from planar diagrams, plus the Skein.jl
# interchange format (signed integer codes).
#
# The walk follows the cycles of the diagram's slot graph (see
# `_slot_graph`): arriving at a slot, the strand passes through the crossing
# along its strand pair (slots 1-3 under, slots 2-4 over) and continues on
# the arc label of the exit slot. Which pair is used decides the over/under
# character of the visit; no local orientation assumption is required.

"""
    GaussVisit(crossing::Int, over::Bool, sign::Int)

One encounter of a traversal with a crossing: `crossing` is the crossing
index, `over` records whether the traversing strand passes OVER, and `sign`
is the +1/-1 crossing sign.
"""
struct GaussVisit
    crossing::Int
    over::Bool
    sign::Int
end

"""
    GaussWord(visits::Vector{GaussVisit}, components::Vector{UnitRange{Int}})

A Gauss word: the sequence of `GaussVisit`s met while walking the link
component cycles, together with the index ranges of each component in
`visits`.
"""
struct GaussWord
    visits::Vector{GaussVisit}
    components::Vector{UnitRange{Int}}
end

"""
    gauss_word(d::PlanarDiagram) -> GaussWord

Walk the diagram and record every crossing encounter. For a knot the word
has `2n` visits, each crossing appearing once as an over-pass and once as
an under-pass.
"""
function gauss_word(d::PlanarDiagram)::GaussWord
    isempty(d.crossings) && return GaussWord(GaussVisit[], UnitRange{Int}[])
    g = _slot_graph(d)
    visited = Set{Int}()
    visits = GaussVisit[]
    ranges = UnitRange{Int}[]
    for slot0 in 1:length(g.label)
        slot0 in visited && continue
        lo = length(visits) + 1
        slot = slot0
        while !(slot in visited)
            push!(visited, slot)
            push!(visited, g.through[slot])
            i = div(slot - 1, 4) + 1
            pos = mod1(slot, 4)
            cr = d.crossings[i]
            push!(visits, GaussVisit(i, pos == 2 || pos == 4, cr.sign))
            slot = g.partner[g.through[slot]]
        end
        push!(ranges, lo:length(visits))
    end
    GaussWord(visits, ranges)
end

"""
    skein_gauss_code(d::PlanarDiagram) -> Vector{Int}

The Skein.jl interchange format for the diagram's Gauss code: a signed
integer sequence where `|k|` identifies the `k`-th crossing in diagram
order and the sign records over (+) versus under (-) at each visit. Every
crossing therefore appears exactly twice, once with each sign — the shape
validated by `Skein.GaussCode`.
"""
function skein_gauss_code(d::PlanarDiagram)::Vector{Int}
    isempty(d.crossings) && return Int[]
    gw = gauss_word(d)
    [v.over ? v.crossing : -v.crossing for v in gw.visits]
end

function Base.show(io::IO, gw::GaussWord)
    code = [v.over ? v.crossing : -v.crossing for v in gw.visits]
    print(io, "GaussWord(", code, ")")
end

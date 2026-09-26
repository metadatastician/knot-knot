# SPDX-License-Identifier: MPL-2.0
# diagram-kernel.jl — minimal working prototype for the diagram kernel
# (src/diagrams.jl + src/gauss.jl).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/diagram-kernel.jl`; it exits
# non-zero unless every pinned assertion holds.
#
# What this proves: the planar-diagram kernel is implementable as pure
# functions over immutable tuples, in the KnotAtlas crossing convention
# shared with KnotTheory.jl:
#
#     X[a, b, c, d] with sign s
#
# the UNDER-strand enters at edge `a` and exits at edge `c`; the OVER-strand
# enters at edge `d` and exits at edge `b`; `s` is +1 for a positive
# (right-hand) crossing and -1 for a negative one. Every edge label occurs
# in exactly two slots, so an "edge" is a strand segment joining two
# crossings.
#
# The pins below are literature/table values, not self-derived: the
# right-hand trefoil 3_1 is all-positive (writhe +3) and the figure-eight
# 4_1 is alternating with two crossings of each sign (writhe 0).

const FAILURES = Ref(0)

function check(label, got, want)
    if got == want
        println("  PASS  $label")
    else
        FAILURES[] += 1
        println("  FAIL  $label: got $got, want $want")
    end
end

# --- the data ---------------------------------------------------------------
# One crossing: (a, b, c, d, sign) in the KnotAtlas convention above.
const Crossing = Tuple{Int,Int,Int,Int,Int}

const TREFOIL = Crossing[
    (1, 4, 2, 5,  1),
    (3, 6, 4, 1,  1),
    (5, 2, 6, 3,  1),
]

const FIGURE_EIGHT = Crossing[
    (4, 2, 5, 1, -1),
    (8, 6, 1, 5, -1),
    (6, 3, 7, 4,  1),
    (2, 7, 3, 8,  1),
]

const UNKNOT = Crossing[]

# Positive Hopf link: two components, two positive crossings. It exists to
# pin component counting — a knot-only test set cannot distinguish "one
# component" from "always one component".
const HOPF = Crossing[
    (1, 2, 3, 4,  1),
    (3, 4, 1, 2,  1),
]

# --- kernel: the operations src/diagrams.jl will expose ---------------------

writhe(d) = sum(c[5] for c in d)

crossing_count(d) = length(d)

mirror(d) = [(a, b, c, dd, -s) for (a, b, c, dd, s) in d]

"""Every edge label of the diagram, sorted."""
function arcs_of(d)
    seen = Int[]
    for c in d
        for edge in (c[1], c[2], c[3], c[4])
            edge in seen || push!(seen, edge)
        end
    end
    sort!(seen)
end

"""
    is_valid_pd(d) -> Bool

A well-formed planar diagram: every edge label occurs in exactly two slots
(no dangling edge, no triple incidence). This is the invariant the kernel
must enforce on construction, because every downstream invariant silently
misbehaves on a malformed code.
"""
function is_valid_pd(d)
    counts = Dict{Int,Int}()
    for c in d
        for edge in (c[1], c[2], c[3], c[4])
            counts[edge] = get(counts, edge, 0) + 1
        end
    end
    all(==(2), values(counts))
end

"""
    slot_table(d) -> (partners, slots)

`slots[edge]` lists the two `(crossing, position)` occurrences of an edge;
`partners[(k, pos)]` is the other slot of the same passage inside crossing
`k` — positions 1<->3 (the under-passage a->c) and 4<->2 (the over-passage
d->b).
"""
function slot_table(d)
    slots = Dict{Int,Vector{Tuple{Int,Int}}}()
    for (k, c) in enumerate(d)
        for pos in 1:4
            edge = c[pos]
            push!(get!(slots, edge, Tuple{Int,Int}[]), (k, pos))
        end
    end
    partners = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    for k in eachindex(d)
        partners[(k, 1)] = (k, 3)
        partners[(k, 3)] = (k, 1)
        partners[(k, 4)] = (k, 2)
        partners[(k, 2)] = (k, 4)
    end
    (partners, slots)
end

"""
    gauss_word(d) -> Vector{Int}

The oriented Gauss word: walk the diagram once, recording `+k` every time
the strand passes OVER crossing `k` and `-k` every time it passes UNDER.
A valid knot diagram yields a word in which every index appears exactly
once positive and once negative.
"""
function gauss_word(d)
    isempty(d) && return Int[]
    partners, slots = slot_table(d)
    start = (1, 1)
    word = Int[]
    slot = start
    while true
        k, pos = slot
        push!(word, pos == 2 || pos == 4 ? k : -k)
        slot = partners[slot]
        occurrences = slots[d[slot[1]][slot[2]]]
        slot = occurrences[1] == slot ? occurrences[2] : occurrences[1]
        slot == start && break
    end
    word
end

"""
    component_count(d) -> Int

Number of link components, counted by exhausting the walk: each component
is one closed walk through the slot graph. Both slots of a passage are
marked, otherwise a single walk would cover only half the slots and every
knot would be miscounted as two components.
"""
function component_count(d)
    isempty(d) && return 0
    partners, slots = slot_table(d)
    seen = Set{Tuple{Int,Int}}()
    components = 0
    for candidate in keys(slots)
        for first in slots[candidate]
            first in seen && continue
            components += 1
            slot = first
            while !(slot in seen)
                push!(seen, slot)
                partner = partners[slot]
                push!(seen, partner)
                occurrences = slots[d[partner[1]][partner[2]]]
                slot = occurrences[1] == partner ? occurrences[2] : occurrences[1]
            end
        end
    end
    components
end

# --- the pins ---------------------------------------------------------------

function main()
    println("diagram-kernel prototype")
    println("  (KnotAtlas PD convention: under a->c, over d->b, sign +-1)")

    println("\n[1] well-formedness and counts")
    check("trefoil is a valid PD", is_valid_pd(TREFOIL), true)
    check("figure-eight is a valid PD", is_valid_pd(FIGURE_EIGHT), true)
    check("a broken PD is rejected", is_valid_pd([(1, 2, 3, 4, 1)]), false)
    check("trefoil crossing count", crossing_count(TREFOIL), 3)
    check("figure-eight crossing count", crossing_count(FIGURE_EIGHT), 4)
    check("trefoil edge count", length(arcs_of(TREFOIL)), 6)
    check("figure-eight edge count", length(arcs_of(FIGURE_EIGHT)), 8)

    println("\n[2] writhe (table values)")
    check("trefoil writhe (all positive)", writhe(TREFOIL), 3)
    check("figure-eight writhe (alternating)", writhe(FIGURE_EIGHT), 0)
    check("mirrored trefoil writhe", writhe(mirror(TREFOIL)), -3)
    check("unknot writhe", writhe(UNKNOT), 0)

    println("\n[3] Gauss words")
    tre = gauss_word(TREFOIL)
    fe = gauss_word(FIGURE_EIGHT)
    check("trefoil Gauss word", tre, [-1, 3, -2, 1, -3, 2])
    check("figure-eight Gauss word", fe, [-1, 2, -3, 4, -2, 1, -4, 3])
    for (name, word) in (("trefoil", tre), ("figure-eight", fe))
        ks = unique(abs.(word))
        balanced = all(k -> count(==(k), word) == 1 && count(==(-k), word) == 1, ks)
        check("$name: each crossing once over, once under", balanced, true)
        check("$name: Gauss length is twice the crossings",
              length(word), 2 * crossing_count(name == "trefoil" ? TREFOIL : FIGURE_EIGHT))
    end
    check("unknot Gauss word is empty", gauss_word(UNKNOT), Int[])

    println("\n[4] components")
    check("trefoil is one component", component_count(TREFOIL), 1)
    check("figure-eight is one component", component_count(FIGURE_EIGHT), 1)
    check("positive Hopf link is two components", component_count(HOPF), 2)
    check("positive Hopf link writhe", writhe(HOPF), 2)
    check("Hopf link is a valid PD", is_valid_pd(HOPF), true)
    check("Hopf link edge count", length(arcs_of(HOPF)), 4)

    println()
    if FAILURES[] == 0
        println("ALL PASS — the diagram kernel prototype is sound.")
    else
        println("$(FAILURES[]) FAILURE(S) — the kernel prototype is wrong.")
        exit(1)
    end
end

main()

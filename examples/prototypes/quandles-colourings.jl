# SPDX-License-Identifier: MPL-2.0
# quandles-colourings.jl — minimal working prototype for the quandle /
# colouring kernel (src/quandles.jl).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/quandles-colourings.jl`; it
# exits non-zero unless every pinned assertion holds.
#
# What this proves: Fox n-colourings are exactly the homomorphisms from the
# fundamental quandle of a diagram into the dihedral quandle R_n = Z/n with
# operation `a #> b = 2b - a`, and the prototype counts them, distinguishes
# the trivial ones, and exhibits the trefoil's rainbow.
#
# The load-bearing detail: a colouring is a function on the ARCS of the
# diagram — the segments that run from one undercrossing to the next — and an
# arc may pass over several crossings. In the PD notation `X[a,b,c,d,s]` the
# over strand runs from `d` to `b`, so `d` and `b` denote the SAME arc and
# must be merged before any colouring is counted. Using the 2c PD labels as
# variables instead of the c arcs over-counts by a factor of n^c (the trefoil
# would report 27 colourings for n = 3 instead of 9).
#
# The colouring rule at a crossing, independent of the sign, is
#     2 * over - under_in = under_out   (mod n)
# with under_in = a, over = d, under_out = c.

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

const PD = Tuple{Int,Int,Int,Int,Int}

const UNKNOT = PD[]

const TREFOIL = PD[
    (1, 4, 2, 5,  1),
    (3, 6, 4, 1,  1),
    (5, 2, 6, 3,  1),
]

const FIGURE_EIGHT = PD[
    (4, 2, 5, 1, -1),
    (8, 6, 1, 5, -1),
    (6, 3, 7, 4,  1),
    (2, 7, 3, 8,  1),
]

const CINQUEFOIL = PD[
    (1, 6, 2, 7,  1),
    (3, 8, 4, 9,  1),
    (5, 10, 6, 1,  1),
    (7, 2, 8, 3,  1),
    (9, 4, 10, 5,  1),
]

# --- the dihedral quandle R_n ------------------------------------------------

"""The dihedral quandle R_n: `a #> b = 2b - a (mod n)`."""
dihedral_op(a, b, n) = mod(2 * b - a, n)

"""Check the three quandle axioms for R_n."""
function quandle_axioms(n)
    elts = 0:(n - 1)
    # a #> a = a
    idem = all(dihedral_op(a, a, n) == a for a in elts)
    # for every a, b there is a unique c with c #> b = a
    right_inv = true
    for a in elts, b in elts
        cs = [c for c in elts if dihedral_op(c, b, n) == a]
        length(cs) == 1 || (right_inv = false)
    end
    # (a #> b) #> c = (a #> c) #> (b #> c)
    self_dist = true
    for a in elts, b in elts, c in elts
        lhs = dihedral_op(dihedral_op(a, b, n), c, n)
        rhs = dihedral_op(dihedral_op(a, c, n), dihedral_op(b, c, n), n)
        lhs == rhs || (self_dist = false)
    end
    (idem, right_inv, self_dist)
end

# --- arcs of the diagram ----------------------------------------------------

mutable struct UnionFind
    parent::Vector{Int}
    rank::Vector{Int}
end

UnionFind(n::Int) = UnionFind(collect(1:n), zeros(Int, n))

function find!(uf::UnionFind, x::Int)
    while uf.parent[x] != x
        uf.parent[x] = uf.parent[uf.parent[x]]
        x = uf.parent[x]
    end
    x
end

function union!(uf::UnionFind, x::Int, y::Int)
    rx = find!(uf, x)
    ry = find!(uf, y)
    rx == ry && return
    if uf.rank[rx] < uf.rank[ry]
        rx, ry = ry, rx
    end
    uf.parent[ry] = rx
    uf.rank[rx] += uf.rank[rx] == uf.rank[ry]
    nothing
end

"""
    arc_map(d) -> (Int, Dict{Int,Int})

The arcs of the diagram: PD labels merged across over-crossings. Returns the
number of arcs and the label -> arc (1-based) map. For a knot diagram with
`c` crossings this is exactly `c`.
"""
function arc_map(d)
    labels = sort(unique(Iterators.flatten(((c[1], c[2], c[3], c[4]) for c in d))))
    idx = Dict(a => i for (i, a) in enumerate(labels))
    uf = UnionFind(length(labels))
    for c in d
        _, b, _, dd, _ = c
        union!(uf, idx[dd], idx[b])          # the over strand: d and b are one arc
    end
    seen = Dict{Int,Int}()
    out = Dict{Int,Int}()
    for a in labels
        r = find!(uf, idx[a])
        if !haskey(seen, r)
            seen[r] = length(seen) + 1
        end
        out[a] = seen[r]
    end
    length(seen), out
end

"""
    colouring_equations(d) -> Vector{Tuple{Int,Int,Int}}

One `(under_in_arc, over_arc, under_out_arc)` triple per crossing: the
constraint is `2 * over - under_in = under_out (mod n)`.
"""
function colouring_equations(d)
    _, am = arc_map(d)
    [(am[c[1]], am[c[4]], am[c[3]]) for c in d]
end

# --- counting colourings ----------------------------------------------------

"""
    count_colourings(d, n) -> Int

Number of Fox n-colourings, i.e. homomorphisms from the fundamental quandle
of the diagram into R_n. Brute force over the arcs: fine for the small
diagrams pinned here, and the honest reference for a smarter search.
"""
function count_colourings(d, n)
    k, _ = arc_map(d)
    k == 0 && return n                          # the unknot has one free arc
    eqs = colouring_equations(d)
    count = 0
    for assign in Iterators.product(ntuple(_ -> 0:(n - 1), k)...)
        ok = true
        for (ui, ov, uo) in eqs
            if mod(2 * assign[ov] - assign[ui], n) != mod(assign[uo], n)
                ok = false
                break
            end
        end
        ok && (count += 1)
    end
    count
end

"""Every colouring, as a tuple of arc colours (1-based arc order)."""
function all_colourings(d, n)
    k, _ = arc_map(d)
    k == 0 && return [()]                       # nothing to colour
    eqs = colouring_equations(d)
    out = []
    for assign in Iterators.product(ntuple(_ -> 0:(n - 1), k)...)
        ok = true
        for (ui, ov, uo) in eqs
            if mod(2 * assign[ov] - assign[ui], n) != mod(assign[uo], n)
                ok = false
                break
            end
        end
        ok && push!(out, assign)
    end
    out
end

"""
    is_n_colourable(d, n) -> Bool

Non-trivially colourable: more than the `n` monochromatic colourings exist.
"""
is_n_colourable(d, n) = count_colourings(d, n) > n

# --- the pins ---------------------------------------------------------------

function main()
    println("quandles-colourings prototype")
    println("  (Fox n-colourings = homomorphisms into the dihedral quandle R_n)")

    println("\n[1] R_n really is a quandle")
    for n in (3, 4, 5)
        idem, right_inv, self_dist = quandle_axioms(n)
        check("R_$n is idempotent", idem, true)
        check("R_$n is right-invertible", right_inv, true)
        check("R_$n is self-distributive", self_dist, true)
    end

    println("\n[2] arc count equals the crossing count")
    check("unknot arcs", arc_map(UNKNOT)[1], 0)
    check("trefoil arcs", arc_map(TREFOIL)[1], 3)
    check("figure-eight arcs", arc_map(FIGURE_EIGHT)[1], 4)
    check("cinquefoil arcs", arc_map(CINQUEFOIL)[1], 5)

    println("\n[3] Fox 3-colourings")
    check("unknot R3", count_colourings(UNKNOT, 3), 3)
    check("trefoil R3", count_colourings(TREFOIL, 3), 9)
    check("figure-eight R3", count_colourings(FIGURE_EIGHT, 3), 3)
    check("cinquefoil R3", count_colourings(CINQUEFOIL, 3), 3)

    println("\n[4] Fox 4-colourings and Fox 5-colourings")
    check("unknot R4", count_colourings(UNKNOT, 4), 4)
    check("trefoil R4", count_colourings(TREFOIL, 4), 4)
    check("figure-eight R4", count_colourings(FIGURE_EIGHT, 4), 4)
    check("cinquefoil R4", count_colourings(CINQUEFOIL, 4), 4)
    check("unknot R5", count_colourings(UNKNOT, 5), 5)
    check("trefoil R5", count_colourings(TREFOIL, 5), 5)
    check("figure-eight R5", count_colourings(FIGURE_EIGHT, 5), 25)
    check("cinquefoil R5", count_colourings(CINQUEFOIL, 5), 25)

    println("\n[5] non-trivial colourability")
    check("unknot is not 3-colourable", is_n_colourable(UNKNOT, 3), false)
    check("trefoil is 3-colourable", is_n_colourable(TREFOIL, 3), true)
    check("figure-eight is not 3-colourable", is_n_colourable(FIGURE_EIGHT, 3), false)
    check("cinquefoil is not 3-colourable", is_n_colourable(CINQUEFOIL, 3), false)
    check("trefoil is not 5-colourable", is_n_colourable(TREFOIL, 5), false)
    check("figure-eight is 5-colourable", is_n_colourable(FIGURE_EIGHT, 5), true)
    check("cinquefoil is 5-colourable", is_n_colourable(CINQUEFOIL, 5), true)

    println("\n[6] the trefoil's 3-colourings are the rainbow")
    # Every Fox 3-colouring of the trefoil is either monochromatic or a
    # permutation of (0, 1, 2): 3 trivial + 6 rainbow = 9.
    trivial = [(c, c, c) for c in 0:2]
    rainbow = [p for p in Iterators.product(0:2, 0:2, 0:2) if length(Set(p)) == 3]
    expected = Set(vcat(trivial, rainbow))
    got = Set(all_colourings(TREFOIL, 3))
    check("trefoil 3-colourings = trivial + rainbow", got, expected)
    check("rainbow has 6 members", length(rainbow), 6)
    # the rainbow really is non-trivial: it uses all three colours
    check("some trefoil 3-colouring uses 3 colours",
          any(c -> length(Set(c)) == 3, collect(got)), true)
end

main()

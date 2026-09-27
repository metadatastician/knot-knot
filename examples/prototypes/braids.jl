# SPDX-License-Identifier: MPL-2.0
# braids.jl — minimal working prototype for the braid kernel (src/braid.jl).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/braids.jl`; it exits non-zero
# unless every pinned assertion holds.
#
# What this proves: a braid word can be closed into a planar diagram in the
# KnotAtlas PD convention, and the closure really is the link it should be —
# checked three independent ways that must agree:
#
#   1. the link components of the closure PD (cycles of the through-strand
#      pairing a<->c, d<->b) equal the cycles of the permutation the braid
#      induces on its strands;
#   2. the determinant, computed exactly from the Fox colouring matrix of the
#      closure (minor of the c x c matrix, one row and one column dropped),
#      matches the literature value of the link;
#   3. the genus from the Seifert circles of the closure matches the
#      (k - 1)/2 formula for the (2,k)-torus knots.
#
# Conventions mirror src/braid.jl: generators are `(i, sign)` with `sigma_i`
# = sign +1 (strand at position i passes OVER the strand at i+1); closure
# connects bottom strand k back to top strand k. A positive generator closes
# to `X[in_(i+1), in_i, out_i, out_(i+1)]`, a negative one to
# `X[in_i, out_i, out_(i+1), in_(i+1)]`.

const FAILURES = Ref(0)

function check(label, got, want)
    if got == want
        println("  PASS  $label")
    else
        FAILURES[] += 1
        println("  FAIL  $label: got $got, want $want")
    end
end

# --- braid words ------------------------------------------------------------

"""A braid word as `(strand_index, sign)` pairs, `+1` for sigma_i and `-1` for its inverse."""
struct Braid
    generators::Vector{Tuple{Int,Int}}
end

"""Parse TANGLE-style braid words: `s1.s1.s1` for positive, `S2` for inverse."""
function parse_braid_word(word::AbstractString)
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

strand_count(b::Braid) = isempty(b.generators) ? 0 : maximum(g[1] for g in b.generators) + 1
writhe(b::Braid) = sum(s for (_, s) in b.generators)

# --- closing the braid ------------------------------------------------------

"""
    braid_closure(b) -> Vector{Tuple{Tuple{Int,Int,Int,Int},Int}}

The planar diagram of the closure, as `(arcs, sign)` pairs. The empty braid
closes to the unknot (no crossings).
"""
function braid_closure(b::Braid)
    isempty(b.generators) && return Tuple{Tuple{Int,Int,Int,Int},Int}[]
    n = strand_count(b)
    current_arc = collect(1:n)
    next_arc = n + 1
    crossings = Tuple{Tuple{Int,Int,Int,Int},Int}[]
    for (i, sgn) in b.generators
        in_i = current_arc[i]
        in_i1 = current_arc[i + 1]
        out_i = next_arc
        out_i1 = next_arc + 1
        next_arc += 2
        if sgn > 0
            push!(crossings, ((in_i1, in_i, out_i, out_i1), 1))
        else
            push!(crossings, ((in_i, out_i, out_i1, in_i1), -1))
        end
        current_arc[i] = out_i
        current_arc[i + 1] = out_i1
    end
    rename = Dict{Int,Int}()
    for k in 1:n
        current_arc[k] != k && (rename[current_arc[k]] = k)
    end
    isempty(rename) && return crossings
    [(Tuple(get(rename, a, a) for a in arcs), s) for (arcs, s) in crossings]
end

# --- the diagram invariants -------------------------------------------------

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

diagram_labels(d) = sort(unique(Iterators.flatten((c[1] for c in d))))

"""
    arc_map(d) -> (Int, Dict{Int,Int})

The arcs of the diagram: PD labels merged across over-crossings, because the
over strand of `X[a,b,c,d]` joins `d` to `b`. For the closure of a knot braid
this is exactly the crossing count.
"""
function arc_map(d)
    labels = diagram_labels(d)
    idx = Dict(a => i for (i, a) in enumerate(labels))
    uf = UnionFind(length(labels))
    for (arcs, _) in d
        _, b, _, dd = arcs
        union!(uf, idx[dd], idx[b])
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

"""Link components: the cycles of the through-strand pairing a<->c, d<->b."""
function link_components(d)
    isempty(d) && return 1                     # the unknot
    labels = diagram_labels(d)
    idx = Dict(a => i for (i, a) in enumerate(labels))
    uf = UnionFind(length(labels))
    for (arcs, _) in d
        a, b, cc, dd = arcs
        union!(uf, idx[a], idx[cc])
        union!(uf, idx[dd], idx[b])
    end
    length(Set(find!(uf, i) for i in eachindex(labels)))
end

"""Seifert circles of the diagram (oriented smoothing), for the genus."""
function seifert_circles(d)
    isempty(d) && return 0
    labels = diagram_labels(d)
    idx = Dict(a => i for (i, a) in enumerate(labels))
    uf = UnionFind(length(labels))
    for (arcs, s) in d
        a, b, cc, dd = arcs
        if s >= 0
            union!(uf, idx[a], idx[dd])
            union!(uf, idx[b], idx[cc])
        else
            union!(uf, idx[a], idx[b])
            union!(uf, idx[cc], idx[dd])
        end
    end
    length(Set(find!(uf, i) for i in eachindex(labels)))
end

genus(d) = isempty(d) ? 0 : div(length(d) - seifert_circles(d) + 1, 2)

"""
    braid_permutation(b) -> Vector{Int}

The permutation of the strands induced by the braid. The closure's components
are the cycles of this permutation.
"""
function braid_permutation(b::Braid)
    n = strand_count(b)
    p = collect(1:n)
    for (i, _) in b.generators
        p[i], p[i + 1] = p[i + 1], p[i]
    end
    p
end

"""Number of cycles of a permutation."""
function cycle_count(p::Vector{Int})
    seen = falses(length(p))
    count = 0
    for k in eachindex(p)
        seen[k] && continue
        count += 1
        j = k
        while !seen[j]
            seen[j] = true
            j = p[j]
        end
    end
    count
end

"""
    determinant(d) -> Int

The determinant of a link diagram, computed exactly from its Fox colouring
matrix: rows are crossings, columns are arcs, `+2` on the over arc and `-1`
on each under arc. Delete the last row and last column and take the absolute
value of the determinant. A one-crossing knot has an empty minor, whose
determinant is 1.
"""
function determinant(d)
    k, am = arc_map(d)
    (k < 2 || length(d) < 2) && return 1
    M = zeros(Int, length(d), k)
    for (r, (arcs, _)) in enumerate(d)
        a, b, cc, dd = arcs
        M[r, am[dd]] += 2
        M[r, am[a]] -= 1
        M[r, am[cc]] -= 1
    end
    A = [Rational{Int}(M[i, j]) for i in 1:(length(d) - 1), j in 1:(k - 1)]
    det = Rational{Int}(1)
    row = 1
    for col in 1:size(A, 2)
        piv = 0
        for i in row:size(A, 1)
            if !iszero(A[i, col])
                piv = i
                break
            end
        end
        piv == 0 && return 0
        if piv != row
            for j in 1:size(A, 2)
                A[row, j], A[piv, j] = A[piv, j], A[row, j]
            end
            det = -det
        end
        det *= A[row, col]
        for i in (row + 1):size(A, 1)
            iszero(A[i, col]) && continue
            f = A[i, col] / A[row, col]
            for j in col:size(A, 2)
                A[i, j] -= f * A[row, j]
            end
        end
        row += 1
    end
    abs(Int(numerator(det)))
end

mirror_word(b::Braid) = Braid([(i, -s) for (i, s) in b.generators])

# --- the pins ---------------------------------------------------------------

function main()
    println("braids prototype")
    println("  (braid words, closure to PD, and three independent agreement checks)")

    println("\n[1] parsing TANGLE-style words")
    check("\"s1.s1.s1\" parses to three sigma_1",
          parse_braid_word("s1.s1.s1").generators, [(1, 1), (1, 1), (1, 1)])
    check("\"S2\" parses to sigma_2^-1",
          parse_braid_word("S2").generators, [(2, -1)])
    check("whitespace and empty tokens are tolerated",
          parse_braid_word("  s1 . . s2  ").generators, [(1, 1), (2, 1)])
    check("the empty word is the empty braid",
          parse_braid_word("").generators, Tuple{Int,Int}[])
    ok = false
    try
        parse_braid_word("x1")
    catch
        ok = true
    end
    check("an invalid generator raises", ok, true)
    ok = false
    try
        parse_braid_word("s0")
    catch
        ok = true
    end
    check("a zero generator index raises", ok, true)

    println("\n[2] closure of the two-strand family (2,k)-torus links")
    check("empty braid closes to the unknot",
          (length(braid_closure(parse_braid_word(""))), link_components(braid_closure(parse_braid_word("")))), (0, 1))
    check("sigma_1 closes to the unknot",
          (length(braid_closure(parse_braid_word("s1"))),
           link_components(braid_closure(parse_braid_word("s1"))),
           determinant(braid_closure(parse_braid_word("s1")))), (1, 1, 1))
    check("sigma_1^2 closes to the Hopf link",
          (length(braid_closure(parse_braid_word("s1.s1"))),
           link_components(braid_closure(parse_braid_word("s1.s1"))),
           determinant(braid_closure(parse_braid_word("s1.s1")))), (2, 2, 2))
    check("sigma_1^3 closes to the trefoil",
          (length(braid_closure(parse_braid_word("s1.s1.s1"))),
           link_components(braid_closure(parse_braid_word("s1.s1.s1"))),
           determinant(braid_closure(parse_braid_word("s1.s1.s1")))), (3, 1, 3))
    check("sigma_1^5 closes to the cinquefoil",
          (length(braid_closure(parse_braid_word("s1.s1.s1.s1.s1"))),
           link_components(braid_closure(parse_braid_word("s1.s1.s1.s1.s1"))),
           determinant(braid_closure(parse_braid_word("s1.s1.s1.s1.s1")))), (5, 1, 5))
    check("sigma_1^-3 closes to the left trefoil",
          (length(braid_closure(parse_braid_word("S1.S1.S1"))),
           link_components(braid_closure(parse_braid_word("S1.S1.S1"))),
           determinant(braid_closure(parse_braid_word("S1.S1.S1")))), (3, 1, 3))

    println("\n[3] a redundant diagram of the same link")
    # sigma_1^2 sigma_1^-1 = sigma_1, so this closes to the unknot with a
    # kink: three crossings, one component, determinant 1.
    red = parse_braid_word("s1.s1.S1")
    check("sigma_1^2 sigma_1^-1 closes to the unknot",
          (length(braid_closure(red)), link_components(braid_closure(red)),
           determinant(braid_closure(red))), (3, 1, 1))
    check("its writhe is +1 (the kink is not cancelled in the word)",
          writhe(red), 1)

    println("\n[4] the three-strand figure-eight")
    fig8 = parse_braid_word("s1.S2.s1.S2")
    d = braid_closure(fig8)
    check("figure-eight closure: 4 crossings, 3 strands, 1 component",
          (length(d), strand_count(fig8), link_components(d)), (4, 3, 1))
    check("figure-eight determinant", determinant(d), 5)
    check("figure-eight braid writhe", writhe(fig8), 0)

    println("\n[5] components agree with the strand permutation")
    for word in ("s1", "s1.s1", "s1.s1.s1", "s1.s1.s1.s1.s1", "S1.S1.S1",
                 "s1.s1.S1", "s1.S2.s1.S2", "s1.s2.s1", "s2.s1.s2",
                 "s1.s2.s1.s2.s1.s2")
        b = parse_braid_word(word)
        check("$word: closure components = permutation cycles",
              link_components(braid_closure(b)), cycle_count(braid_permutation(b)))
    end

    println("\n[6] arcs equal crossings for knot closures")
    for word in ("s1", "s1.s1.s1", "s1.s1.s1.s1.s1", "s1.S2.s1.S2")
        b = parse_braid_word(word)
        check("$word: arcs = crossings", arc_map(braid_closure(b))[1], length(braid_closure(b)))
    end

    println("\n[7] genus of the (2,k)-torus knots")
    check("trefoil genus", genus(braid_closure(parse_braid_word("s1.s1.s1"))), 1)
    check("cinquefoil genus", genus(braid_closure(parse_braid_word("s1.s1.s1.s1.s1"))), 2)
    check("figure-eight genus", genus(braid_closure(parse_braid_word("s1.S2.s1.S2"))), 1)

    println("\n[8] conjugation and mirroring")
    # The half-twists sigma_1 sigma_2 sigma_1 and sigma_2 sigma_1 sigma_2 are
    # conjugate, and both close to the Hopf link - a different braid with the
    # same closure as sigma_1^2.
    for word in ("s1.s2.s1", "s2.s1.s2")
        b = parse_braid_word(word)
        check("$word closes to the Hopf link",
              (link_components(braid_closure(b)), determinant(braid_closure(b))), (2, 2))
    end
    # mirroring the word mirrors the closure: same determinant
    for word in ("s1.s1.s1", "s1.s1.s1.s1.s1", "s1.S2.s1.S2")
        b = parse_braid_word(word)
        check("mirror of $word has the same determinant",
              determinant(braid_closure(mirror_word(b))),
              determinant(braid_closure(b)))
    end

    println()
    if FAILURES[] == 0
        println("ALL PASS — the braid-kernel prototype is sound.")
    else
        println("$(FAILURES[]) FAILURE(S) — the braid-kernel prototype is wrong.")
        exit(1)
    end
end

main()

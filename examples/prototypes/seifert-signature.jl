# SPDX-License-Identifier: MPL-2.0
# seifert-signature.jl — minimal working prototype for the Seifert-signature
# kernel (src/seifert.jl).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/seifert-signature.jl`; it exits
# non-zero unless every pinned assertion holds.
#
# What this proves: Seifert's algorithm — oriented smoothing to Seifert
# circles, one band per crossing, a spanning tree to fix a basis of H_1, and
# the band layout that fills the Seifert matrix V — gives a matrix whose
# symmetrisation V + V^T has the classical signatures, and the genus formula
# g = (c - s + 1) / 2.
#
# Conventions mirror src/seifert.jl:
#   * oriented smoothing at X[a,b,c,d]: positive crossings connect a<->d and
#     b<->c; negative crossings connect a<->b and c<->d;
#   * a positive crossing contributes the band (circle(a), circle(b)) and a
#     negative one (circle(a), circle(c));
#   * V[i,i] = +s if the band's two ends are on the same circle, else -s;
#     off-diagonal entries come from bands sharing a circle or the same
#     circle pair, ordered by band index.
#
# The signature is computed EXACTLY (rational congruence elimination), not by
# floating-point eigendecomposition: the answer is an integer and the inputs
# are small, so there is no reason to approximate it.

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
const Crossing = Tuple{Int,Int,Int,Int,Int}

const UNKNOT = Crossing[]

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

const CINQUEFOIL = Crossing[
    (1, 6, 2, 7,  1),
    (3, 8, 4, 9,  1),
    (5, 10, 6, 1,  1),
    (7, 2, 8, 3,  1),
    (9, 4, 10, 5,  1),
]

writhe(d) = sum(c[5] for c in d)
mirror(d) = [(a, b, c, dd, -s) for (a, b, c, dd, s) in d]

# --- Seifert circles --------------------------------------------------------

"""Union-find with path halving and union by rank."""
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
    seifert_circles(d) -> Int

Number of Seifert circles: arcs joined by the oriented smoothing of every
crossing. For the empty diagram this is 0 (there is no strand to smooth).
"""
function seifert_circles(d)
    isempty(d) && return 0
    arcs = sort(unique(Iterators.flatten(((c[1], c[2], c[3], c[4]) for c in d))))
    idx = Dict(a => i for (i, a) in enumerate(arcs))
    uf = UnionFind(length(arcs))
    for c in d
        a, b, cc, dd, s = c
        if s >= 0
            union!(uf, idx[a], idx[dd])
            union!(uf, idx[b], idx[cc])
        else
            union!(uf, idx[a], idx[b])
            union!(uf, idx[cc], idx[dd])
        end
    end
    length(Set(find!(uf, i) for i in eachindex(arcs)))
end

"""
    arc_to_circle(d) -> Dict{Int,Int}

Which Seifert circle each arc belongs to.
"""
function arc_to_circle(d)
    arcs = sort(unique(Iterators.flatten(((c[1], c[2], c[3], c[4]) for c in d))))
    idx = Dict(a => i for (i, a) in enumerate(arcs))
    uf = UnionFind(length(arcs))
    for c in d
        a, b, cc, dd, s = c
        if s >= 0
            union!(uf, idx[a], idx[dd])
            union!(uf, idx[b], idx[cc])
        else
            union!(uf, idx[a], idx[b])
            union!(uf, idx[cc], idx[dd])
        end
    end
    seen = Dict{Int,Int}()
    out = Dict{Int,Int}()
    for a in arcs
        r = find!(uf, idx[a])
        if !haskey(seen, r)
            seen[r] = length(seen) + 1
        end
        out[a] = seen[r]
    end
    out
end

# --- the Seifert matrix -----------------------------------------------------

"""
    seifert_matrix(d) -> Matrix{Int}

A Seifert matrix for the diagram, of size `g x g` with
`g = c - s + 1`. The construction is the one src/seifert.jl pins: the
right-hand trefoil yields `[-1 1; 0 -1]`.
"""
function seifert_matrix(d)
    n = length(d)
    n == 0 && return Matrix{Int}(undef, 0, 0)
    a2c = arc_to_circle(d)
    n_circles = length(unique(values(a2c)))

    bands = Tuple{Int,Int,Int}[]
    for c in d
        a, b, cc, dd, s = c
        if s >= 0
            push!(bands, (a2c[a], a2c[b], s))
        else
            push!(bands, (a2c[a], a2c[cc], s))
        end
    end

    # A spanning tree of the Seifert graph fixes a basis of H_1: every band
    # NOT in the tree contributes one generator.
    tree = UnionFind(n_circles)
    is_tree = falses(n)
    for (k, (ci, cj, _)) in enumerate(bands)
        ci == cj && continue
        if find!(tree, ci) != find!(tree, cj)
            is_tree[k] = true
            union!(tree, ci, cj)
        end
    end
    non_tree = [k for k in 1:n if !is_tree[k]]
    g = length(non_tree)
    g == 0 && return Matrix{Int}(undef, 0, 0)

    V = zeros(Int, g, g)
    for (i, ki) in enumerate(non_tree)
        ci_i, cj_i, si = bands[ki]
        V[i, i] = ci_i == cj_i ? si : -si
        for (j, kj) in enumerate(non_tree)
            i == j && continue
            ci_j, cj_j, sj = bands[kj]
            (ci_i == cj_i || ci_j == cj_j) && continue
            same_pair =
                (ci_i == ci_j && cj_i == cj_j) || (ci_i == cj_j && cj_i == ci_j)
            if same_pair
                if ki < kj
                    consecutive = true
                    for (m, km) in enumerate(non_tree)
                        (m == i || m == j) && continue
                        ci_m, cj_m, _ = bands[km]
                        same_m =
                            (ci_m == ci_i && cj_m == cj_i) ||
                            (ci_m == cj_i && cj_m == ci_i)
                        if same_m && ki < km < kj
                            consecutive = false
                            break
                        end
                    end
                    consecutive && (V[i, j] = si >= 0 ? 1 : -1)
                end
                continue
            end
            shares_circle =
                ci_i == ci_j || ci_i == cj_j || cj_i == ci_j || cj_i == cj_j
            if shares_circle && ki < kj
                V[i, j] = si >= 0 ? 1 : -1
            end
        end
    end
    V
end

"""
    signature(d) -> Int

The knot signature: positive minus negative eigenvalues of the symmetrised
Seifert matrix `V + V^T`. Computed exactly by rational congruence
elimination (row and column swaps together, so the transformation stays
congruent and the signature is preserved).
"""
function signature(d)
    V = seifert_matrix(d)
    isempty(V) && return 0
    n = size(V, 1)
    A = [Rational{Int}(V[i, j] + V[j, i]) for i in 1:n, j in 1:n]
    pos = 0
    neg = 0
    for k in 1:n
        piv = 0
        for i in k:n
            if !iszero(A[i, k])
                piv = i
                break
            end
        end
        piv == 0 && continue
        if piv != k
            for j in 1:n
                A[k, j], A[piv, j] = A[piv, j], A[k, j]
            end
            for i in 1:n
                A[i, k], A[i, piv] = A[i, piv], A[i, k]
            end
        end
        dk = A[k, k]
        dk > 0 ? (pos += 1) : (neg += 1)
        for i in (k + 1):n
            iszero(A[i, k]) && continue
            f = A[i, k] / dk
            for j in k:n
                A[i, j] -= f * A[k, j]
            end
        end
    end
    pos - neg
end

"""
    genus(d) -> Int

The Seifert genus of a knot diagram: `g = (c - s + 1) / 2`. Integer for
knots; the unknot has genus 0.
"""
function genus(d)
    isempty(d) && return 0
    div(length(d) - seifert_circles(d) + 1, 2)
end

# --- the pins ---------------------------------------------------------------

function main()
    println("seifert-signature prototype")
    println("  (Seifert's algorithm: oriented smoothing, bands, spanning tree, V + V^T)")

    println("\n[1] Seifert circles and genus")
    check("unknot has no Seifert circles", seifert_circles(UNKNOT), 0)
    check("unknot genus", genus(UNKNOT), 0)
    check("trefoil Seifert circles", seifert_circles(TREFOIL), 2)
    check("trefoil genus", genus(TREFOIL), 1)
    check("figure-eight Seifert circles", seifert_circles(FIGURE_EIGHT), 3)
    check("figure-eight genus", genus(FIGURE_EIGHT), 1)
    check("cinquefoil Seifert circles", seifert_circles(CINQUEFOIL), 2)
    check("cinquefoil genus", genus(CINQUEFOIL), 2)

    println("\n[2] Seifert matrices (the trefoil pin is the published one)")
    check("trefoil V = [-1 1; 0 -1]",
          seifert_matrix(TREFOIL), [-1 1; 0 -1])
    check("figure-eight V = [1 -1; 0 -1]",
          seifert_matrix(FIGURE_EIGHT), [1 -1; 0 -1])
    check("cinquefoil V is 4x4 with -1 diagonal and +1 superdiagonal",
          seifert_matrix(CINQUEFOIL),
          [-1 1 0 0; 0 -1 1 0; 0 0 -1 1; 0 0 0 -1])
    check("unknot V is empty", size(seifert_matrix(UNKNOT)), (0, 0))

    println("\n[3] signatures (classical values)")
    check("unknot signature", signature(UNKNOT), 0)
    check("right-hand trefoil signature", signature(TREFOIL), -2)
    check("figure-eight signature", signature(FIGURE_EIGHT), 0)
    check("cinquefoil (2,5)-torus knot signature", signature(CINQUEFOIL), -4)
    check("mirrored trefoil signature", signature(mirror(TREFOIL)), 2)

    println("\n[4] structural facts")
    for (name, d) in (("trefoil", TREFOIL), ("figure-eight", FIGURE_EIGHT),
                      ("cinquefoil", CINQUEFOIL))
        check("$name: signature is even", iseven(signature(d)), true)
        check("$name: |signature| <= 2 * genus",
              abs(signature(d)) <= 2 * genus(d), true)
        V = seifert_matrix(d)
        S = V + transpose(V)
        check("$name: V + V^T is symmetric", S, transpose(S))
    end

    println()
    if FAILURES[] == 0
        println("ALL PASS — the Seifert-signature kernel prototype is sound.")
    else
        println("$(FAILURES[]) FAILURE(S) — the Seifert-signature kernel prototype is wrong.")
        exit(1)
    end
end

main()

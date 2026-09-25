# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Seifert surface data: Seifert circles, the band-based Seifert matrix, the
# signature, and the genus.
#
# Oriented smoothing at a crossing X[a,b,c,d] follows the KnotAtlas
# convention: positive crossings connect a<->d and b<->c; negative
# crossings connect a<->b and c<->d. The Seifert graph has one vertex per
# circle and one edge (band) per crossing; a spanning tree fixes a basis of
# H_1 and the band layout fills the Seifert matrix V. The construction is
# pinned by the right-hand trefoil, which yields V = [-1 1; 0 -1].

"""
    seifert_circles(d::PlanarDiagram) -> Int

Number of Seifert circles obtained by oriented smoothing of every crossing.
"""
function seifert_circles(d::PlanarDiagram)::Int
    isempty(d.crossings) && return 0
    all_arcs = arcs_of(d)
    idx = Dict{Int,Int}(a => i for (i, a) in enumerate(all_arcs))
    parent = collect(1:length(all_arcs))
    rnk = zeros(Int, length(all_arcs))
    for c in d.crossings
        a, b, cc, dd = c.arcs
        if c.sign >= 0
            _uf_union!(parent, rnk, idx[a], idx[dd])
            _uf_union!(parent, rnk, idx[b], idx[cc])
        else
            _uf_union!(parent, rnk, idx[a], idx[b])
            _uf_union!(parent, rnk, idx[cc], idx[dd])
        end
    end
    length(Set(_uf_find!(parent, i) for i in eachindex(all_arcs)))
end

"""
    seifert_matrix(d::PlanarDiagram) -> Matrix{Int}

A Seifert matrix ``V`` for the diagram, of size ``g \\times g`` with
``g = c - s + 1`` (crossings minus Seifert circles plus one). Right-hand
trefoil yields ``[-1\\ 1; 0\\ -1]``.
"""
function seifert_matrix(d::PlanarDiagram)::Matrix{Int}
    n = length(d.crossings)
    if n == 0
        return Matrix{Int}(undef, 0, 0)
    end

    # Circle labels per arc via the same smoothing union-find.
    all_arcs = arcs_of(d)
    idx = Dict{Int,Int}(a => i for (i, a) in enumerate(all_arcs))
    parent = collect(1:length(all_arcs))
    rnk = zeros(Int, length(all_arcs))
    for c in d.crossings
        a, b, cc, dd = c.arcs
        if c.sign >= 0
            _uf_union!(parent, rnk, idx[a], idx[dd])
            _uf_union!(parent, rnk, idx[b], idx[cc])
        else
            _uf_union!(parent, rnk, idx[a], idx[b])
            _uf_union!(parent, rnk, idx[cc], idx[dd])
        end
    end
    root_to_circle = Dict{Int,Int}()
    arc_to_circle = Dict{Int,Int}()
    n_circles = 0
    for a in all_arcs
        r = _uf_find!(parent, idx[a])
        if !haskey(root_to_circle, r)
            n_circles += 1
            root_to_circle[r] = n_circles
        end
        arc_to_circle[a] = root_to_circle[r]
    end

    bands = Tuple{Int,Int,Int}[]
    for c in d.crossings
        a, b, cc, dd = c.arcs
        if c.sign >= 0
            push!(bands, (arc_to_circle[a], arc_to_circle[b], c.sign))
        else
            push!(bands, (arc_to_circle[a], arc_to_circle[cc], c.sign))
        end
    end

    tree_parent = collect(1:n_circles)
    tree_rank = zeros(Int, n_circles)
    is_tree = falses(n)
    for (k, (ci, cj, _)) in enumerate(bands)
        ci == cj && continue
        if _uf_find!(tree_parent, ci) != _uf_find!(tree_parent, cj)
            is_tree[k] = true
            _uf_union!(tree_parent, tree_rank, ci, cj)
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
            same_pair = (ci_i == ci_j && cj_i == cj_j) || (ci_i == cj_j && cj_i == ci_j)
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
    signature(d::PlanarDiagram) -> Int

The knot signature ``\\sigma(K)``: the signature (positive minus negative
eigenvalue count) of the symmetrised Seifert matrix ``V + V^T``.
Right-hand trefoil: ``-2``; figure-eight: ``0``.
"""
function signature(d::PlanarDiagram)::Int
    V = seifert_matrix(d)
    isempty(V) && return 0
    S = Symmetric(Float64.(V + transpose(V)))
    ev = eigvals(S)
    count(x -> x > 1e-10, ev) - count(x -> x < -1e-10, ev)
end

"""
    genus(d::PlanarDiagram) -> Int

The Seifert genus of a knot diagram: ``g = (c - s + 1) / 2`` where ``c``
is the crossing number and ``s`` the Seifert circle count. For knots this
is an integer; the unknot has genus 0.
"""
function genus(d::PlanarDiagram)::Int
    isempty(d.crossings) && return 0
    s = seifert_circles(d)
    c = length(d.crossings)
    div(c - s + 1, 2)
end

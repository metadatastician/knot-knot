# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The diagram kernel: KnotAtlas-convention planar-diagram (PD) codes.
#
# Convention (shared with KnotTheory.jl, which pins it against KnotAtlas):
# a crossing is recorded as `(a, b, c, d, sign)` with the four incident arcs
# listed counter-clockwise. The UNDER-strand enters at arc `a` and exits at
# arc `c`; the OVER-strand enters at arc `d` and exits at arc `b`. The sign
# is +1 for a positive (right-hand) crossing and -1 for a negative one.

"""
    Crossing(arcs::NTuple{4,Int}, sign::Int)

One crossing of a planar diagram in the KnotAtlas convention described in
the module docstring. `sign` is +1 (positive) or -1 (negative).
"""
struct Crossing
    arcs::NTuple{4,Int}
    sign::Int
end

"""
    PlanarDiagram(crossings::Vector{Crossing}, components::Vector{Vector{Int}})

An oriented link diagram given by its crossings and, optionally, a grouping
of arc labels by link component. For knots `components` may be empty.
"""
struct PlanarDiagram
    crossings::Vector{Crossing}
    components::Vector{Vector{Int}}
end

"""
    pdcode(entries::Vector{NTuple{5,Int}}; components = Vector{Vector{Int}}())

Build a `PlanarDiagram` from raw `(a, b, c, d, sign)` crossing tuples.

# Example — right-hand trefoil (KnotAtlas X[1,4,2,5], X[3,6,4,1], X[5,2,6,3])
```julia
pd = pdcode([(1, 4, 2, 5, 1), (3, 6, 4, 1, 1), (5, 2, 6, 3, 1)])
```
"""
function pdcode(
    entries::Vector{NTuple{5,Int}};
    components::Vector{Vector{Int}} = Vector{Vector{Int}}(),
)::PlanarDiagram
    PlanarDiagram([Crossing((e[1], e[2], e[3], e[4]), e[5]) for e in entries], components)
end

"""
    pd(x) -> Vector{NTuple{5,Int}}

Extract the `(a, b, c, d, sign)` tuples of a `PlanarDiagram`.
"""
pd(d::PlanarDiagram)::Vector{NTuple{5,Int}} =
    [(c.arcs[1], c.arcs[2], c.arcs[3], c.arcs[4], c.sign) for c in d.crossings]

"""
    crossing_count(d::PlanarDiagram) -> Int

Number of crossings of the diagram.
"""
crossing_count(d::PlanarDiagram)::Int = length(d.crossings)

"""
    writhe(d::PlanarDiagram) -> Int

Sum of the crossing signs. The writhe is a diagram invariant, not a knot
invariant (it changes under Reidemeister I moves).
"""
writhe(d::PlanarDiagram)::Int = isempty(d.crossings) ? 0 : sum(c.sign for c in d.crossings)

"""
    mirror(d::PlanarDiagram) -> PlanarDiagram

Mirror image: flip every crossing sign (leave the arc layout untouched).
For knot invariants the mirror sends ``V(t) \\mapsto V(t^{-1})`` and
``\\sigma \\mapsto -\\sigma``.
"""
function mirror(d::PlanarDiagram)::PlanarDiagram
    PlanarDiagram([Crossing(c.arcs, -c.sign) for c in d.crossings], d.components)
end

"""
    arc_positions(d::PlanarDiagram) -> Dict{Int,Vector{Int}}

Map each arc label to the crossing slots where it appears. Slot numbering is
`4 * (crossing_index - 1) + position`, position in 1:4, so every slot is a
globally unique integer usable by union-find based state sums.
"""
function arc_positions(d::PlanarDiagram)::Dict{Int,Vector{Int}}
    pos = Dict{Int,Vector{Int}}()
    for (i, c) in enumerate(d.crossings)
        for (slot, arc) in enumerate(c.arcs)
            push!(get!(pos, arc, Int[]), 4 * (i - 1) + slot)
        end
    end
    pos
end

# --- small union-find used by several kernels ------------------------------

function _uf_find!(parent::Vector{Int}, x::Int)::Int
    while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end
    x
end

function _uf_union!(parent::Vector{Int}, rnk::Vector{Int}, x::Int, y::Int)
    rx, ry = _uf_find!(parent, x), _uf_find!(parent, y)
    rx == ry && return
    if rnk[rx] < rnk[ry]
        parent[rx] = ry
    elseif rnk[rx] > rnk[ry]
        parent[ry] = rx
    else
        parent[ry] = rx
        rnk[rx] += 1
    end
    nothing
end

"""
    arcs_of(d::PlanarDiagram) -> Vector{Int}

Sorted list of every arc label used by the diagram.
"""
function arcs_of(d::PlanarDiagram)::Vector{Int}
    s = Set{Int}()
    for c in d.crossings
        for a in c.arcs
            push!(s, a)
        end
    end
    sort(collect(s))
end

"""
    wirt_generators(d::PlanarDiagram) -> (gen_map::Dict{Int,Int}, n_gens::Int)

Wirtinger generators: at each crossing the two over-strand arcs `b` and `d`
belong to one generator (they are one continuous strand of the knot). Returns
a map from arc label to 1-based generator index and the number of generators.
"""
function wirt_generators(d::PlanarDiagram)::Tuple{Dict{Int,Int},Int}
    isempty(d.crossings) && return (Dict{Int,Int}(), 0)
    all_arcs = arcs_of(d)
    idx = Dict{Int,Int}(a => i for (i, a) in enumerate(all_arcs))
    parent = collect(1:length(all_arcs))
    rnk = zeros(Int, length(all_arcs))
    for c in d.crossings
        _uf_union!(parent, rnk, idx[c.arcs[4]], idx[c.arcs[2]])
    end
    root_to_gen = Dict{Int,Int}()
    gen_map = Dict{Int,Int}()
    n = 0
    for a in all_arcs
        r = _uf_find!(parent, idx[a])
        if !haskey(root_to_gen, r)
            n += 1
            root_to_gen[r] = n
        end
        gen_map[a] = root_to_gen[r]
    end
    (gen_map, n)
end

# Slot graph of a diagram. Slots are numbered 4 * (crossing - 1) + position.
# Two kinds of edges: THROUGH edges join the two slots of one strand pair at
# a crossing (slots 1-3 carry the under-strand, slots 2-4 the over-strand),
# and ARC edges join the two occurrences of one arc label. The graph is
# 2-regular, so its cycles are exactly the link components — no orientation
# assumption is needed (PD codes do not fix strand directions locally).
function _slot_graph(d::PlanarDiagram)
    n_slots = 4 * length(d.crossings)
    label = zeros(Int, n_slots)
    occ = Dict{Int,Vector{Int}}()
    for (i, c) in enumerate(d.crossings)
        for (s, a) in enumerate(c.arcs)
            slot = 4 * (i - 1) + s
            label[slot] = a
            push!(get!(occ, a, Int[]), slot)
        end
    end
    for (a, v) in occ
        length(v) == 2 || error("arc $a appears $(length(v)) times: not a diagram")
    end
    partner = zeros(Int, n_slots)
    for v in values(occ)
        partner[v[1]] = v[2]
        partner[v[2]] = v[1]
    end
    through = zeros(Int, n_slots)
    for slot in 1:n_slots
        pos = mod1(slot, 4)
        through[slot] = slot + (pos == 1 ? 2 : pos == 2 ? 2 : pos == 3 ? -2 : -2)
    end
    (label = label, partner = partner, through = through)
end

"""
    traverse_components(d::PlanarDiagram) -> Vector{Vector{Int}}

Walk every cycle of the diagram's slot graph, returning one arc sequence per
link component. At a slot the strand passes THROUGH the crossing to the
other slot of its pair (slots 1-3: under, slots 2-4: over) and then travels
along that arc label to its other occurrence. Throws if any arc label does
not occur exactly twice.
"""
function traverse_components(d::PlanarDiagram)::Vector{Vector{Int}}
    isempty(d.crossings) && return Vector{Int}[]
    g = _slot_graph(d)
    visited = Set{Int}()
    comps = Vector{Int}[]
    for slot0 in 1:length(g.label)
        slot0 in visited && continue
        arcs = Int[]
        slot = slot0
        while !(slot in visited)
            push!(visited, slot)
            push!(visited, g.through[slot])
            push!(arcs, g.label[slot])
            slot = g.partner[g.through[slot]]
        end
        push!(comps, arcs)
    end
    comps
end

"""
    slot_orientation(d::PlanarDiagram) -> Vector{Int}

A consistent strand orientation of the diagram's slots, as a 0/1 vector
indexed by global slot number: 1 = the strand ENTERS the crossing at this
slot, 0 = it exits. Opposite slots of one crossing always carry opposite
values, as do the two occurrences of one arc label, so the orientation is
a 2-colouring of the slot constraint graph — bipartite for every valid
link diagram (the constructor throws otherwise). The colouring is fixed by
choosing the first slot of each connected component as an entry; reversing
a component's orientation flips all of its values.
"""
function slot_orientation(d::PlanarDiagram)::Vector{Int}
    n_slots = 4 * length(d.crossings)
    isempty(d.crossings) && return Int[]
    adj = [Int[] for _ in 1:n_slots]
    occ = Dict{Int,Vector{Int}}()
    for (i, c) in enumerate(d.crossings)
        for (s, a) in enumerate(c.arcs)
            push!(get!(occ, a, Int[]), 4 * (i - 1) + s)
        end
    end
    for v in values(occ)
        length(v) == 2 || error("arc appears $(length(v)) times: not a diagram")
        push!(adj[v[1]], v[2])
        push!(adj[v[2]], v[1])
    end
    for i in 0:(length(d.crossings) - 1)
        b = 4 * i
        push!(adj[b + 1], b + 3)
        push!(adj[b + 3], b + 1)
        push!(adj[b + 2], b + 4)
        push!(adj[b + 4], b + 2)
    end
    color = fill(-1, n_slots)
    for s0 in 1:n_slots
        color[s0] != -1 && continue
        color[s0] = 1
        stack = [s0]
        while !isempty(stack)
            u = pop!(stack)
            for w in adj[u]
                if color[w] == -1
                    color[w] = 1 - color[u]
                    push!(stack, w)
                elseif color[w] == color[u]
                    error("diagram is not orientable")
                end
            end
        end
    end
    color
end

"""
    linking_number(d::PlanarDiagram) -> Int

The pairwise linking number of a 2-component diagram: half the signed sum
of the crossings where the two components meet. Throws for diagrams that do
not have exactly two components.
"""
function linking_number(d::PlanarDiagram)::Int
    comps = traverse_components(d)
    length(comps) == 2 || error("linking_number requires a 2-component diagram")
    comp_of = Dict{Int,Int}()
    for (k, arcs) in enumerate(comps), a in arcs
        comp_of[a] = k
    end
    s = 0
    for c in d.crossings
        comp_of[c.arcs[1]] != comp_of[c.arcs[2]] || continue
        s += c.sign
    end
    div(s, 2)
end

function Base.show(io::IO, d::PlanarDiagram)
    print(io, "PlanarDiagram(crossings=$(length(d.crossings)), writhe=$(writhe(d)))")
end

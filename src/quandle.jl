# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Finite quandles and the fundamental quandle of a diagram.
#
# A quandle is a set Q with an operation |> satisfying:
#   (Q1) idempotence        x |> x = x
#   (Q2) right-invertibility for each y, x -> x |> y is a bijection
#   (Q3) self-distributivity (x |> y) |> z = (x |> z) |> (y |> z)
#
# At a crossing of a diagram the Wirtinger arcs x (under, entering),
# z (under, exiting) and y (over) satisfy z = x |> y at a POSITIVE crossing
# and z = x |>^-1 y at a NEGATIVE one. The number of homomorphisms from the
# fundamental quandle into a finite quandle Q is the Q-colouring count, a
# knot invariant.

"""
    Quandle(table::Matrix{Int})

A finite quandle on the elements `1:n`, given by its Cayley table:
`table[x, y]` is ``x \\triangleright y``. Construct through
`dihedral_quandle` / `alexander_quandle` or supply a table checked by
`is_quandle`.
"""
struct Quandle
    table::Matrix{Int}
end

Base.length(q::Quandle)::Int = size(q.table, 1)

operate(q::Quandle, x::Int, y::Int)::Int = q.table[x, y]

"""
    is_quandle(q::Quandle) -> Bool

Check the three quandle axioms against the Cayley table.
"""
function is_quandle(q::Quandle)::Bool
    n = length(q)
    size(q.table) == (n, n) || return false
    for x in 1:n
        q.table[x, x] == x || return false
    end
    for y in 1:n
        sort([q.table[x, y] for x in 1:n]) == collect(1:n) || return false
    end
    for x in 1:n, y in 1:n, z in 1:n
        q.table[q.table[x, y], z] == q.table[q.table[x, z], y] || return false
    end
    true
end

"""
    dihedral_quandle(p::Int) -> Quandle

The dihedral (Alexander) quandle ``R_p`` on ``\\mathbb{Z}/p``:
``x \\triangleright y = 2y - x \\pmod p``. ``R_3``-colourability is
classical tricolourability.
"""
function dihedral_quandle(p::Int)::Quandle
    p >= 1 || error("dihedral_quandle: p must be positive")
    table = Matrix{Int}(undef, p, p)
    for x in 1:p, y in 1:p
        table[x, y] = mod(2(y - 1) - (x - 1), p) + 1
    end
    Quandle(table)
end

"""
    alexander_quandle(n::Int, t::Int) -> Quandle

The Alexander quandle on ``\\mathbb{Z}/n``: ``x \\triangleright y =
t x + (1 - t) y \\pmod n``. A quandle for every unit `t` mod `n`.
"""
function alexander_quandle(n::Int, t::Int)::Quandle
    n >= 1 || error("alexander_quandle: n must be positive")
    table = Matrix{Int}(undef, n, n)
    for x in 1:n, y in 1:n
        table[x, y] = mod(t * (x - 1) + (1 - t) * (y - 1), n) + 1
    end
    Quandle(table)
end

"""
    fundamental_relations(d::PlanarDiagram) -> Vector{Tuple{Int,Int,Int,Int}}

The fundamental quandle presentation of the diagram as
`(x, y, z, sign)` tuples: under-generator `x` enters the crossing,
over-generator `y` crosses it, under-generator `z` exits, and `sign` is the
crossing sign (+1: ``z = x \\triangleright y``; -1: ``z = x \\triangleright^{-1} y``).
Generators are the Wirtinger arcs.
"""
function fundamental_relations(d::PlanarDiagram)::Vector{Tuple{Int,Int,Int,Int}}
    isempty(d.crossings) && return Tuple{Int,Int,Int,Int}[]
    gen_map, _ = wirt_generators(d)
    orient = slot_orientation(d)
    rels = Tuple{Int,Int,Int,Int}[]
    for (i, c) in enumerate(d.crossings)
        a, b, cc, dd = c.arcs
        base = 4 * (i - 1)
        x, z = orient[base + 1] == 1 ? (gen_map[a], gen_map[cc]) : (gen_map[cc], gen_map[a])
        push!(rels, (x, gen_map[dd], z, c.sign))
    end
    rels
end

# Solve z = x |> y for x, using column-bijectivity of the Cayley table.
function _inv_operate(q::Quandle, z::Int, y::Int)::Int
    for x in 1:length(q)
        q.table[x, y] == z && return x
    end
    error("quandle table violates right-invertibility")
end

"""
    coloring_count(d::PlanarDiagram, q::Quandle) -> Int

The number of quandle colourings: homomorphisms from the fundamental quandle
of `d` into `q`. Pinned values: unknot 3, trefoil 9, figure-eight 3 for
``R_3``.
"""
function coloring_count(d::PlanarDiagram, q::Quandle)::Int
    _, n_gens = wirt_generators(d)
    n_gens == 0 && return length(q)  # unknot: any colour
    rels = fundamental_relations(d)
    count = 0
    assignment = zeros(Int, n_gens)

    function ok_at_depth(depth::Int)::Bool
        for (x, y, z, sgn) in rels
            assignment[x] == 0 && continue
            assignment[y] == 0 && continue
            assignment[z] == 0 && continue
            if sgn >= 0
                q.table[assignment[x], assignment[y]] == assignment[z] || return false
            else
                _inv_operate(q, assignment[z], assignment[y]) == assignment[x] || return false
            end
        end
        true
    end

    function backtrack(depth::Int)
        if depth > n_gens
            count += 1
            return
        end
        for v in 1:length(q)
            assignment[depth] = v
            if ok_at_depth(depth)
                backtrack(depth + 1)
            end
        end
        assignment[depth] = 0
    end

    backtrack(1)
    count
end

function Base.show(io::IO, q::Quandle)
    print(io, "Quandle(order=$(length(q)), axioms=$(is_quandle(q)))")
end

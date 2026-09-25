# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Standard knot table (Rolfsen names 0_1 .. 7_7).
#
# Every knot 5_2 .. 7_7 is pinned directly from the KnotAtlas PD skeletons
# (katlas.org, the PD4Knots data used by KnotTheory`) with an explicitly
# verified sign assignment for each crossing. The signs were pinned by
# exhaustive search: the unique assignment whose Kauffman-bracket Jones
# polynomial AND Fox-calculus Alexander polynomial agree exactly with the
# values published on the knot's KnotAtlas page (cross-checked against the
# closure of the Atlas braid word for each knot). The older DT-based route
# (KnotTheory.jl's own DT table plus a signed Gauss reconstruction) does
# NOT reproduce the standard knots beyond 5_1 — its 6-crossing and higher
# entries reconstruct knots with wrong determinants — so no table knot here
# is constructed from a DT code. `from_dt` remains available as a utility,
# but the table itself never goes through it.
#
# Provenance per knot: PD skeleton from the KnotAtlas `PD[K]` listing;
# signs pinned against the Atlas Alexander/Jones polynomials; determinants,
# signatures and genera below are the Atlas "Determinant and Signature" and
# "3-genus" values, all reproduced by this package's own implementations in
# the test suite.

"""
    KnotEntry(name::String, pd::PlanarDiagram, dt::Vector{Int}, description::String)

One row of the standard knot table. `dt` records the KnotAtlas
Dowker-Thistlethwaite code for reference only; the diagram `pd` is built
from the pinned signed PD, never from the DT code.
"""
struct KnotEntry
    name::String
    pd::PlanarDiagram
    dt::Vector{Int}
    description::String
end

const _TREFOIL_PD = PlanarDiagram(
    [
        Crossing((1, 4, 2, 5), 1),
        Crossing((3, 6, 4, 1), 1),
        Crossing((5, 2, 6, 3), 1),
    ],
    Vector{Vector{Int}}(),
)

const _FIGURE_EIGHT_PD = PlanarDiagram(
    [
        Crossing((4, 2, 5, 1), -1),
        Crossing((8, 6, 1, 5), -1),
        Crossing((6, 3, 7, 4), 1),
        Crossing((2, 7, 3, 8), 1),
    ],
    Vector{Vector{Int}}(),
)

const _CINQUEFOIL_PD = PlanarDiagram(
    [
        Crossing((1, 6, 2, 7), 1),
        Crossing((3, 8, 4, 9), 1),
        Crossing((5, 10, 6, 1), 1),
        Crossing((7, 2, 8, 3), 1),
        Crossing((9, 4, 10, 5), 1),
    ],
    Vector{Vector{Int}}(),
)

# --- pinned KnotAtlas PDs for 5_2 .. 7_7 ----------------------------------
# Each block: the Atlas skeleton (X[a,b,c,d] per crossing, slots listed
# anticlockwise starting at the entering under-arc) together with the
# crossing signs pinned against the Atlas polynomial invariants.

const _5_2_PD = PlanarDiagram(
    [
        Crossing((1, 4, 2, 5), 1),
        Crossing((3, 8, 4, 9), 1),
        Crossing((5, 10, 6, 1), 1),
        Crossing((9, 6, 10, 7), 1),
        Crossing((7, 2, 8, 3), 1),
    ],
    Vector{Vector{Int}}(),
)

const _6_1_PD = PlanarDiagram(
    [
        Crossing((1, 4, 2, 5), 1),
        Crossing((7, 10, 8, 11), 1),
        Crossing((3, 9, 4, 8), -1),
        Crossing((9, 3, 10, 2), -1),
        Crossing((5, 12, 6, 1), 1),
        Crossing((11, 6, 12, 7), 1),
    ],
    Vector{Vector{Int}}(),
)

const _6_2_PD = PlanarDiagram(
    [
        Crossing((1, 4, 2, 5), 1),
        Crossing((5, 10, 6, 11), 1),
        Crossing((3, 9, 4, 8), 1),
        Crossing((9, 3, 10, 2), 1),
        Crossing((7, 12, 8, 1), -1),
        Crossing((11, 6, 12, 7), -1),
    ],
    Vector{Vector{Int}}(),
)

const _6_3_PD = PlanarDiagram(
    [
        Crossing((4, 2, 5, 1), 1),
        Crossing((8, 4, 9, 3), 1),
        Crossing((12, 9, 1, 10), 1),
        Crossing((10, 5, 11, 6), -1),
        Crossing((6, 11, 7, 12), -1),
        Crossing((2, 8, 3, 7), -1),
    ],
    Vector{Vector{Int}}(),
)

const _7_1_PD = PlanarDiagram(
    [
        Crossing((1, 8, 2, 9), 1),
        Crossing((3, 10, 4, 11), 1),
        Crossing((5, 12, 6, 13), 1),
        Crossing((7, 14, 8, 1), 1),
        Crossing((9, 2, 10, 3), 1),
        Crossing((11, 4, 12, 5), 1),
        Crossing((13, 6, 14, 7), 1),
    ],
    Vector{Vector{Int}}(),
)

const _7_2_PD = PlanarDiagram(
    [
        Crossing((1, 4, 2, 5), 1),
        Crossing((3, 10, 4, 11), 1),
        Crossing((5, 14, 6, 1), 1),
        Crossing((7, 12, 8, 13), 1),
        Crossing((11, 8, 12, 9), 1),
        Crossing((13, 6, 14, 7), 1),
        Crossing((9, 2, 10, 3), 1),
    ],
    Vector{Vector{Int}}(),
)

const _7_3_PD = PlanarDiagram(
    [
        Crossing((6, 2, 7, 1), -1),
        Crossing((10, 4, 11, 3), -1),
        Crossing((14, 8, 1, 7), -1),
        Crossing((8, 14, 9, 13), -1),
        Crossing((12, 6, 13, 5), -1),
        Crossing((2, 10, 3, 9), -1),
        Crossing((4, 12, 5, 11), -1),
    ],
    Vector{Vector{Int}}(),
)

const _7_4_PD = PlanarDiagram(
    [
        Crossing((6, 2, 7, 1), -1),
        Crossing((12, 6, 13, 5), -1),
        Crossing((14, 8, 1, 7), -1),
        Crossing((8, 14, 9, 13), -1),
        Crossing((2, 12, 3, 11), -1),
        Crossing((10, 4, 11, 3), -1),
        Crossing((4, 10, 5, 9), -1),
    ],
    Vector{Vector{Int}}(),
)

const _7_5_PD = PlanarDiagram(
    [
        Crossing((1, 4, 2, 5), 1),
        Crossing((3, 10, 4, 11), 1),
        Crossing((5, 12, 6, 13), 1),
        Crossing((7, 14, 8, 1), 1),
        Crossing((13, 6, 14, 7), 1),
        Crossing((11, 8, 12, 9), 1),
        Crossing((9, 2, 10, 3), 1),
    ],
    Vector{Vector{Int}}(),
)

const _7_6_PD = PlanarDiagram(
    [
        Crossing((1, 4, 2, 5), 1),
        Crossing((3, 8, 4, 9), 1),
        Crossing((5, 12, 6, 13), 1),
        Crossing((9, 1, 10, 14), 1),
        Crossing((13, 11, 14, 10), 1),
        Crossing((11, 6, 12, 7), 1),
        Crossing((7, 2, 8, 3), 1),
    ],
    Vector{Vector{Int}}(),
)

const _7_7_PD = PlanarDiagram(
    [
        Crossing((1, 4, 2, 5), 1),
        Crossing((5, 10, 6, 11), 1),
        Crossing((3, 9, 4, 8), -1),
        Crossing((9, 3, 10, 2), -1),
        Crossing((11, 14, 12, 1), 1),
        Crossing((7, 13, 8, 12), -1),
        Crossing((13, 7, 14, 6), -1),
    ],
    Vector{Vector{Int}}(),
)

# --- DT utilities (kept as an API; the table itself does not use them) ----

# DT -> signed Gauss code: entry i pairs odd position 2i-1 with even
# position |code[i]|; the sign of code[i] decides which occurrence is the
# over-pass (+ => odd position over).
function _dt_to_signed_gauss(code::Vector{Int})::Vector{Int}
    n = length(code)
    n == 0 && return Int[]
    gauss = zeros(Int, 2n)
    for i in 1:n
        odd_pos = 2i - 1
        even_pos = abs(code[i])
        1 <= even_pos <= 2n || error("invalid DT entry $(code[i])")
        iseven(even_pos) || error("DT entry $(code[i]) must pair with an even position")
        gauss[odd_pos] != 0 && error("odd position $odd_pos used twice")
        gauss[even_pos] != 0 && error("even position $even_pos used twice")
        if code[i] >= 0
            gauss[odd_pos] = i
            gauss[even_pos] = -i
        else
            gauss[odd_pos] = -i
            gauss[even_pos] = i
        end
    end
    gauss
end

# Signed Gauss -> PD in the KnotAtlas convention. The odd-position occurrence
# fixes the crossing sign (its sign), and the over/under roles fix the slot
# layout (under in/out at slots 1/3, over in/out at slots 4/2).
function _signed_gauss_to_pd(code::Vector{Int})::PlanarDiagram
    isempty(code) && return PlanarDiagram(Crossing[], Vector{Vector{Int}}())
    n = length(code) ÷ 2
    pos = Dict{Int,Vector{Int}}()
    for (p, t) in enumerate(code)
        push!(get!(pos, abs(t), Int[]), p)
    end
    entries = Vector{NTuple{5,Int}}()
    for i in 1:n
        p1, p2 = pos[i]
        t1, t2 = code[p1], code[p2]
        over_pos, under_pos = if t1 > 0 && t2 < 0
            (p1, p2)
        elseif t2 > 0 && t1 < 0
            (p2, p1)
        else
            (min(p1, p2), max(p1, p2))
        end
        odd_pos = if isodd(p1)
            p1
        elseif isodd(p2)
            p2
        else
            over_pos
        end
        crossing_sign = code[odd_pos] >= 0 ? 1 : -1
        under_in = mod1(under_pos - 1, 2n)
        over_in = mod1(over_pos - 1, 2n)
        push!(entries, (under_in, over_pos, under_pos, over_in, crossing_sign))
    end
    pdcode(entries)
end

"""
    from_dt(code::Vector{Int}) -> PlanarDiagram

Build a planar diagram from a Dowker-Thistlethwaite code via the signed
Gauss-code route. Canonical diagrams are used for 3_1, 4_1 and 5_1 where
hand-verified KnotAtlas PD codes are pinned.

!!! warning
    For knots with more than five crossings, DT reconstruction fixes the
    knot only up to a crossing-sign pattern that the code alone does not
    determine; the reconstruction can land on a different knot with the
    same shadow. The standard table therefore pins every knot 5_2 .. 7_7
    from signed KnotAtlas PDs instead of going through this function.
"""
function from_dt(code::Vector{Int})::PlanarDiagram
    code == Int[] && return PlanarDiagram(Crossing[], Vector{Vector{Int}}())
    code == [4, 6, 2] && return _TREFOIL_PD
    code == [4, 6, 8, 2] && return _FIGURE_EIGHT_PD
    code == [6, 8, 10, 2, 4] && return _CINQUEFOIL_PD
    _signed_gauss_to_pd(_dt_to_signed_gauss(code))
end

# name => (pinned PD, reference DT code from the KnotAtlas, description)
const _KNOT_TABLE_DATA = [
    ("0_1", PlanarDiagram(Crossing[], Vector{Vector{Int}}()), Int[], "Unknot"),
    ("3_1", _TREFOIL_PD, [4, 6, 2], "Right-hand trefoil, the (2,3)-torus knot"),
    ("4_1", _FIGURE_EIGHT_PD, [4, 6, 8, 2], "The figure-eight knot"),
    ("5_1", _CINQUEFOIL_PD, [6, 8, 10, 2, 4], "Cinquefoil, the (2,5)-torus knot"),
    ("5_2", _5_2_PD, [4, 8, 10, 2, 6], "Three-twist knot"),
    ("6_1", _6_1_PD, [4, 8, 12, 10, 2, 6], "Stevedore knot"),
    ("6_2", _6_2_PD, [4, 8, 10, 12, 2, 6], "Miller Institute knot"),
    ("6_3", _6_3_PD, [4, 8, 10, 2, 12, 6], "Six-crossing alternating knot"),
    ("7_1", _7_1_PD, [4, 6, 8, 10, 12, 14, 2], "The (2,7)-torus knot"),
    ("7_2", _7_2_PD, [4, 10, 14, 12, 2, 8, 6], "Seven-crossing twist knot"),
    ("7_3", _7_3_PD, [6, 10, 12, 14, 2, 4, 8], "Seven-crossing alternating knot"),
    ("7_4", _7_4_PD, [6, 10, 12, 14, 4, 2, 8], "Seven-crossing alternating knot"),
    ("7_5", _7_5_PD, [4, 10, 12, 14, 2, 8, 6], "Seven-crossing alternating knot"),
    ("7_6", _7_6_PD, [4, 8, 12, 2, 14, 6, 10], "Seven-crossing alternating knot"),
    ("7_7", _7_7_PD, [4, 8, 10, 12, 2, 14, 6], "Seven-crossing alternating knot"),
]

const _KNOT_TABLE = Dict{String,KnotEntry}()

function _knot_table()::Dict{String,KnotEntry}
    if isempty(_KNOT_TABLE)
        for (name, pd, dt, desc) in _KNOT_TABLE_DATA
            _KNOT_TABLE[name] = KnotEntry(name, pd, dt, desc)
        end
    end
    _KNOT_TABLE
end

"""
    knot_names() -> Vector{String}

Names of every knot in the standard table, sorted by crossing number then
Rolfsen index.
"""
function knot_names()::Vector{String}
    _knot_table()
    sort(collect(keys(_KNOT_TABLE)), by = name -> (length(name), name))
end

"""
    standard_knot(name::AbstractString) -> KnotEntry

Fetch a standard knot by Rolfsen name (`"3_1"`, `"4_1"`, ...; `"0_1"` for
the unknot). Throws `ArgumentError` for unknown names.
"""
function standard_knot(name::AbstractString)::KnotEntry
    t = _knot_table()
    haskey(t, String(name)) || throw(ArgumentError("unknown knot: $name"))
    t[String(name)]
end

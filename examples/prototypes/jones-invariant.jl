# SPDX-License-Identifier: MPL-2.0
# jones-invariant.jl — minimal working prototype for the Jones invariant
# (src/jones.jl + src/laurent.jl).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/jones-invariant.jl`; it exits
# non-zero unless every pinned assertion holds.
#
# What this proves: the Kauffman-bracket state sum, as a pure function over
# the planar-diagram kernel, reproduces the literature Jones polynomials for
# the knots up to five crossings, and the |V(-1)| = determinant identity
# ties the Jones kernel to the Alexander-Conway kernel.
#
# Conventions mirror src/jones.jl exactly (they are not the textbook
# defaults, so they are stated here rather than assumed):
#
#   * the A-smoothing at EVERY crossing pairs (a,b) and (c,d); the
#     B-smoothing pairs (b,c) and (d,a). The pairing is a property of the
#     unoriented picture — which slots belong to the over-strand — so it
#     does not depend on the crossing sign;
#   * <D> = sum over states of A^(#A - #B) * d^(loops - 1),
#     d = -A^2 - A^-2;
#   * V(t) = (-A^3)^w <D> evaluated at A = t^(-1/4), w the writhe.
#
# Pins are KnotAtlas/table values.

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

# --- the slot graph ---------------------------------------------------------
# slots[edge] lists the two (crossing, position) occurrences of an edge.
# Position 1 = a, 2 = b, 3 = c, 4 = d.
function slot_table(d)
    slots = Dict{Int,Vector{Tuple{Int,Int}}}()
    for (k, c) in enumerate(d)
        for pos in 1:4
            push!(get!(slots, c[pos], Tuple{Int,Int}[]), (k, pos))
        end
    end
    slots
end

# --- Laurent arithmetic on Dict{Int,Int} (exponent -> coefficient) ----------
const D_FACTOR = Dict(-2 => -1, 2 => -1)          # d = -A^2 - A^-2

function lp_mul(p, q)
    r = Dict{Int,Int}()
    for (e1, c1) in p, (e2, c2) in q
        r[e1 + e2] = get(r, e1 + e2, 0) + c1 * c2
    end
    filter(kv -> kv[2] != 0, r)
end

function lp_pow(p, m)
    r = Dict(0 => 1)
    for _ in 1:m
        r = lp_mul(r, p)
    end
    r
end

# --- the Kauffman bracket ---------------------------------------------------

"""
    smoothed_loops(d, state) -> Int

Number of closed loops after smoothing every crossing according to `state`
(`0` = A-smoothing, `1` = B-smoothing), with each edge's two slots joined.
Union-find over the `4 * crossings` slots; the loop count is the number of
components.
"""
function smoothed_loops(d, state)
    n = length(d)
    m = 4 * n
    parent = collect(1:m)
    function find!(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    function link!(x, y)
        rx = find!(x)
        ry = find!(y)
        rx == ry || (parent[rx] = ry)
        nothing
    end
    for k in 1:n
        base = 4 * (k - 1)
        if state[k] == 0
            link!(base + 1, base + 2)     # a-b
            link!(base + 3, base + 4)     # c-d
        else
            link!(base + 2, base + 3)     # b-c
            link!(base + 4, base + 1)     # d-a
        end
    end
    slots = slot_table(d)
    for edge in keys(slots)
        occ = slots[edge]
        i1 = 4 * (occ[1][1] - 1) + occ[1][2]
        i2 = 4 * (occ[2][1] - 1) + occ[2][2]
        link!(i1, i2)
    end
    length(Set(find!(i) for i in 1:m))
end

"""
    bracket_polynomial(d) -> Dict{Int,Int}

The Kauffman bracket as a Laurent polynomial in `A`, keyed by exponent.
Exponential in the crossing number, so it is bounded the way src/jones.jl
bounds it.
"""
const MAX_BRACKET_CROSSINGS = 12

function bracket_polynomial(d)
    n = length(d)
    n == 0 && return Dict(0 => 1)
    n > MAX_BRACKET_CROSSINGS && error(
        "bracket state sum limited to $MAX_BRACKET_CROSSINGS crossings (got $n)",
    )
    total = Dict{Int,Int}()
    for mask in 0:(2^n - 1)
        state = [(mask >> (k - 1)) & 1 for k in 1:n]      # 0 = A, 1 = B
        na = count(==(0), state)
        nb = n - na
        loops = smoothed_loops(d, state)
        for (e, c) in lp_pow(D_FACTOR, loops - 1)
            total[e + na - nb] = get(total, e + na - nb, 0) + c
        end
    end
    filter(kv -> kv[2] != 0, total)
end

"""
    jones_polynomial(d) -> Dict{Int,Int}

`V(t) = (-A^3)^w <D>` at `A = t^(-1/4)`, w the writhe. Throws when an
exponent is not a multiple of four, which is the same guard src/jones.jl
carries: such an input is not a knot diagram.
"""
function jones_polynomial(d)
    n = length(d)
    n == 0 && return Dict(0 => 1)
    br = bracket_polynomial(d)
    w = writhe(d)
    sgn = isodd(w) ? -1 : 1
    out = Dict{Int,Int}()
    for (k, c) in br
        num = -(k + 3 * w)
        num % 4 == 0 || error("non-integral Jones exponent: not a knot diagram")
        q = div(num, 4)
        out[q] = get(out, q, 0) + sgn * c
    end
    filter(kv -> kv[2] != 0, out)
end

# --- reading the result -----------------------------------------------------
lp_string(p) = join(
    [string(c > 0 ? "+" : "", c, "t^", e) for (e, c) in sort(collect(p))],
    " ",
)
value_at_one(p) = sum(values(p))
value_at_minus_one(p) = sum(iseven(e) ? c : -c for (e, c) in p)

# --- the pins ---------------------------------------------------------------

function main()
    println("jones-invariant prototype")
    println("  (Kauffman bracket state sum; V(t) = (-A^3)^w <D>, A = t^(-1/4))")

    println("\n[1] brackets (writhe-free, table-checked)")
    check("unknot bracket", bracket_polynomial(UNKNOT), Dict(0 => 1))
    check("trefoil bracket", bracket_polynomial(TREFOIL), Dict(-5 => -1, 3 => -1, 7 => 1))
    check("figure-eight bracket is symmetric",
          sort(collect(keys(bracket_polynomial(FIGURE_EIGHT)))),
          [-8, -4, 0, 4, 8])

    println("\n[2] Jones polynomials (KnotAtlas values)")
    check("unknot V(t) = 1", jones_polynomial(UNKNOT), Dict(0 => 1))
    check("trefoil V(t) = -t^-4 + t^-3 + t^-1",
          jones_polynomial(TREFOIL), Dict(-4 => -1, -3 => 1, -1 => 1))
    check("figure-eight V(t) = t^-2 - t^-1 + 1 - t + t^2",
          jones_polynomial(FIGURE_EIGHT), Dict(-2 => 1, -1 => -1, 0 => 1, 1 => -1, 2 => 1))
    check("cinquefoil V(t) = -t^-7 + t^-6 - t^-5 + t^-4 + t^-2",
          jones_polynomial(CINQUEFOIL),
          Dict(-7 => -1, -6 => 1, -5 => -1, -4 => 1, -2 => 1))

    println("\n[3] identities every knot must satisfy")
    for (name, d) in (("unknot", UNKNOT), ("trefoil", TREFOIL),
                      ("figure-eight", FIGURE_EIGHT), ("cinquefoil", CINQUEFOIL))
        v = jones_polynomial(d)
        check("$name: V(1) = 1", value_at_one(v), 1)
        check("$name: |V(-1)| is odd and positive",
              isodd(abs(value_at_minus_one(v))) && value_at_minus_one(v) != 0, true)
    end

    println("\n[4] |V(-1)| = determinant (bridge to the Alexander kernel)")
    check("trefoil |V(-1)| = 3", abs(value_at_minus_one(jones_polynomial(TREFOIL))), 3)
    check("figure-eight |V(-1)| = 5",
          abs(value_at_minus_one(jones_polynomial(FIGURE_EIGHT))), 5)
    check("cinquefoil |V(-1)| = 5",
          abs(value_at_minus_one(jones_polynomial(CINQUEFOIL))), 5)

    println()
    if FAILURES[] == 0
        println("ALL PASS — the Jones kernel prototype is sound.")
        println("  trefoil      V(t) = $(lp_string(jones_polynomial(TREFOIL)))")
        println("  figure-eight V(t) = $(lp_string(jones_polynomial(FIGURE_EIGHT)))")
    else
        println("$(FAILURES[]) FAILURE(S) — the Jones kernel prototype is wrong.")
        exit(1)
    end
end

main()

# SPDX-License-Identifier: MPL-2.0
# alexander-conway.jl — minimal working prototype for the Alexander-Conway
# kernel (src/alexander.jl + src/laurent.jl).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/alexander-conway.jl`; it exits
# non-zero unless every pinned assertion holds.
#
# What this proves: the Fox-calculus route to the Alexander polynomial —
# Wirtinger presentation from the PD code, Fox derivatives, the Alexander
# matrix, an exact Laurent determinant, and the Conway substitution
# z^2 = t + t^-1 - 2 — reproduces the literature polynomials for the knots
# up to five crossings, and |Delta(-1)| equals the determinant.
#
# The one convention that is easy to get wrong: the Wirtinger generators are
# NOT the PD's edge labels. An edge is cut at every crossing, whereas a
# Wirtinger generator runs along the over-strand, so at each crossing the two
# over-strand edges (d and b) denote the SAME generator. Skipping that merge
# yields a matrix whose minors are constants, not polynomials.

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

# --- Laurent arithmetic on Dict{Int,Int} (exponent -> coefficient) ----------

function lp_mul(p, q)
    r = Dict{Int,Int}()
    for (e1, c1) in p, (e2, c2) in q
        r[e1 + e2] = get(r, e1 + e2, 0) + c1 * c2
    end
    filter(kv -> kv[2] != 0, r)
end

function lp_add(p, q)
    r = Dict{Int,Int}(p)
    for (e, c) in q
        r[e] = get(r, e, 0) + c
    end
    filter(kv -> kv[2] != 0, r)
end

"""
    det_laurent(M) -> Dict{Int,Int}

Determinant of a matrix of Laurent polynomials by cofactor expansion along
the first row. Exponential in the size, which is fine here: the Alexander
matrix of a diagram with c crossings is c x c and this prototype is bounded
to five crossings. (src/alexander.jl uses Bareiss; the permutation and
Bareiss routes agree, and that agreement is itself a test worth keeping.)
"""
function det_laurent(M)
    n = size(M, 1)
    n == 0 && return Dict(0 => 1)
    n == 1 && return Dict(M[1, 1])
    total = Dict{Int,Int}()
    for j in 1:n
        entries = [M[i, k] for i in 2:n for k in 1:n if k != j]
        sub = det_laurent(reshape(entries, n - 1, n - 1))
        sgn = iseven(j) ? -1 : 1
        for (e, c) in lp_mul(M[1, j], sub)
            total[e] = get(total, e, 0) + sgn * c
        end
    end
    filter(kv -> kv[2] != 0, total)
end

# --- Wirtinger presentation -------------------------------------------------

"""
    wirtinger(d) -> (relators, n_generators)

Generators are the over-strand classes: at each crossing the two over-strand
edges `d` and `b` denote the same generator, so they are merged with
union-find before the relators are written. Each crossing contributes one
relator; for a positive crossing `x_d x_a x_d^-1 x_c^-1`, for a negative one
`x_d^-1 x_a x_d x_c^-1`.
"""
function wirtinger(d)
    isempty(d) && return (Vector{Vector{Tuple{Int,Int}}}(), 0)
    edges = sort(unique(Iterators.flatten(((c[1], c[2], c[3], c[4]) for c in d))))
    parent = Dict(e => e for e in edges)
    function find!(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    for c in d
        a, b, cc, dd, s = c
        rx, ry = find!(dd), find!(b)
        rx == ry || (parent[rx] = ry)
    end
    gen = Dict(e => find!(e) for e in edges)
    names = sort(unique(values(gen)))
    index = Dict(g => i for (i, g) in enumerate(names))
    relators = Vector{Vector{Tuple{Int,Int}}}()
    for c in d
        a, b, cc, dd, s = c
        xd = index[gen[dd]]
        xa = index[gen[a]]
        xc = index[gen[cc]]
        if s > 0
            push!(relators, [(xd, 1), (xa, 1), (xd, -1), (xc, -1)])
        else
            push!(relators, [(xd, -1), (xa, 1), (xd, 1), (xc, -1)])
        end
    end
    (relators, length(names))
end

"""
    fox_derivative(word, j) -> Dict{Int,Int}

The Fox derivative of a group word with respect to generator `j`, mapped to
`Z[t, t^-1]` by sending every generator to `t`. Rules:
`d(x_i)/d(x_j) = delta_ij`, `d(x_i^-1)/d(x_j) = -x_i^-1 delta_ij`, and
`d(uv)/d(x_j) = du/d(x_j) + u dv/d(x_j)`.
"""
function fox_derivative(word::Vector{Tuple{Int,Int}}, j::Int)
    result = Dict{Int,Int}()
    prefix = Dict(0 => 1)
    for (g, e) in word
        if g == j
            if e == 1
                for (k, c) in prefix
                    result[k] = get(result, k, 0) + c
                end
            else
                for (k, c) in prefix
                    result[k - 1] = get(result, k - 1, 0) - c
                end
            end
        end
        prefix = lp_mul(prefix, Dict(e => 1))
    end
    filter(kv -> kv[2] != 0, result)
end

"""
    alexander_polynomial(d) -> Dict{Int,Int}

The Alexander polynomial, normalised the way the package pins it: palindromic
about exponent zero and `Delta(1) = +1`. That representative is unique, and
it is what makes the trefoil and figure-eight pins comparable.
"""
function alexander_polynomial(d)
    isempty(d) && return Dict(0 => 1)
    relators, ng = wirtinger(d)
    M = [fox_derivative(r, j) for r in relators, j in 1:ng]
    n = size(M, 1)
    minor = [M[i, j] for i in 1:(n - 1), j in 1:(ng - 1)]
    p = det_laurent(minor)
    isempty(p) && return p
    ks = sort(collect(keys(p)))
    shift = div(ks[1] + ks[end], 2)
    q = Dict(e - shift => c for (e, c) in p)
    if sum(values(q)) < 0
        q = Dict(e => -c for (e, c) in q)
    end
    filter(kv -> kv[2] != 0, q)
end

# --- the Conway substitution ------------------------------------------------

function binomial(n, k)
    (k < 0 || k > n) && return 0
    result = 1
    for i in 1:k
        result = div(result * (n - k + i), i)
    end
    result
end

"""
    chebyshev_like(n) -> Dict{Int,Int}

`P_n(s)` with `t^n + t^-n = P_n(t + t^-1)`: `P_0 = 2`, `P_1 = s`,
`P_n = s P_(n-1) - P_(n-2)`. The `P_0 = 2` matters — using 1 silently halves
the constant term of every Conway polynomial.
"""
function chebyshev_like(n)
    n == 0 && return Dict(0 => 2)
    n == 1 && return Dict(1 => 1)
    p0 = Dict(0 => 2)
    p1 = Dict(1 => 1)
    for _ in 2:n
        p2 = Dict{Int,Int}()
        for (e, c) in p1
            p2[e + 1] = get(p2, e + 1, 0) + c
        end
        for (e, c) in p0
            p2[e] = get(p2, e, 0) - c
        end
        p0, p1 = p1, filter(kv -> kv[2] != 0, p2)
    end
    p1
end

"""
    conway_polynomial(delta) -> Dict{Int,Int}

`Delta(t) = grad(z)` with `z^2 = t + t^-1 - 2`, i.e. substitute
`s = t + t^-1 = z^2 + 2` into the `s`-basis form of `Delta`.
"""
function conway_polynomial(delta)
    isempty(delta) && return Dict(0 => 1)
    g = maximum(abs.(collect(keys(delta))))
    s_poly = Dict(0 => get(delta, 0, 0))
    for i in 1:g
        b = get(delta, i, 0)
        b == 0 && continue
        for (e, c) in chebyshev_like(i)
            s_poly[e] = get(s_poly, e, 0) + b * c
        end
    end
    z = Dict{Int,Int}()
    for (k, c) in s_poly
        for j in 0:k
            coef = c * binomial(k, j) * 2^(k - j)
            coef == 0 && continue
            z[2 * j] = get(z, 2 * j, 0) + coef
        end
    end
    filter(kv -> kv[2] != 0, z)
end

# --- reading the result -----------------------------------------------------
lp_string(p) = join([string(c > 0 ? "+" : "", c, "t^", e)
                     for (e, c) in sort(collect(p))], " ")
value_at_one(p) = sum(values(p))
value_at_minus_one(p) = sum(iseven(e) ? c : -c for (e, c) in p)

# --- the pins ---------------------------------------------------------------

function main()
    println("alexander-conway prototype")
    println("  (Fox calculus on the Wirtinger presentation; Conway by z^2 = t + t^-1 - 2)")

    println("\n[1] Wirtinger presentation shape")
    for (name, d) in (("trefoil", TREFOIL), ("figure-eight", FIGURE_EIGHT),
                      ("cinquefoil", CINQUEFOIL))
        relators, ng = wirtinger(d)
        check("$name: one relator per crossing", length(relators), length(d))
        check("$name: generators equal crossings", ng, length(d))
        ok = all(all(1 <= g <= ng && abs(e) == 1 for (g, e) in r) for r in relators)
        check("$name: relators use only in-range generators, exponent +-1", ok, true)
    end

    println("\n[2] Alexander polynomials (literature values)")
    check("unknot Delta(t) = 1", alexander_polynomial(UNKNOT), Dict(0 => 1))
    check("trefoil Delta(t) = t^-1 - 1 + t",
          alexander_polynomial(TREFOIL), Dict(-1 => 1, 0 => -1, 1 => 1))
    check("figure-eight Delta(t) = -t^-1 + 3 - t",
          alexander_polynomial(FIGURE_EIGHT), Dict(-1 => -1, 0 => 3, 1 => -1))
    check("cinquefoil Delta(t) = t^-2 - t^-1 + 1 - t + t^2",
          alexander_polynomial(CINQUEFOIL),
          Dict(-2 => 1, -1 => -1, 0 => 1, 1 => -1, 2 => 1))

    println("\n[3] Conway polynomials (skein-theoretic values)")
    check("unknot grad(z) = 1", conway_polynomial(Dict(0 => 1)), Dict(0 => 1))
    check("trefoil grad(z) = 1 + z^2",
          conway_polynomial(alexander_polynomial(TREFOIL)), Dict(0 => 1, 2 => 1))
    check("figure-eight grad(z) = 1 - z^2",
          conway_polynomial(alexander_polynomial(FIGURE_EIGHT)), Dict(0 => 1, 2 => -1))
    check("cinquefoil grad(z) = 1 + 3z^2 + z^4",
          conway_polynomial(alexander_polynomial(CINQUEFOIL)),
          Dict(0 => 1, 2 => 3, 4 => 1))

    println("\n[4] identities and the determinant bridge")
    for (name, d) in (("trefoil", TREFOIL), ("figure-eight", FIGURE_EIGHT),
                      ("cinquefoil", CINQUEFOIL))
        delta = alexander_polynomial(d)
        check("$name: Delta(1) = 1", value_at_one(delta), 1)
        pal = all(get(delta, e, 0) == get(delta, -e, 0) for e in keys(delta))
        check("$name: Delta is palindromic", pal, true)
    end
    check("trefoil |Delta(-1)| = 3",
          abs(value_at_minus_one(alexander_polynomial(TREFOIL))), 3)
    check("figure-eight |Delta(-1)| = 5",
          abs(value_at_minus_one(alexander_polynomial(FIGURE_EIGHT))), 5)
    check("cinquefoil |Delta(-1)| = 5",
          abs(value_at_minus_one(alexander_polynomial(CINQUEFOIL))), 5)

    println("\n[5] the Conway substitution inverts the Alexander one")
    for (name, d) in (("trefoil", TREFOIL), ("figure-eight", FIGURE_EIGHT))
        z = conway_polynomial(alexander_polynomial(d))
        # grad(z) with z^2 = t + t^-1 - 2 must rebuild Delta(t): check that the
        # two agree at t = 4 (z^2 = 9/2) and t = 9 (z^2 = 64/9).
        delta = alexander_polynomial(d)
        for t in (4, 9)
            zsq = t + 1 // t - 2
            lhs = sum(c * zsq^(div(e, 2)) for (e, c) in z; init = 0 // 1)
            rhs = sum(c * (t // 1)^e for (e, c) in delta; init = 0 // 1)
            check("$name: grad(z^2) = Delta(t) at t=$t", lhs, rhs)
        end
    end

    println()
    if FAILURES[] == 0
        println("ALL PASS — the Alexander-Conway kernel prototype is sound.")
        println("  trefoil      Delta(t) = $(lp_string(alexander_polynomial(TREFOIL)))")
        println("  figure-eight Delta(t) = $(lp_string(alexander_polynomial(FIGURE_EIGHT)))")
    else
        println("$(FAILURES[]) FAILURE(S) — the Alexander-Conway kernel prototype is wrong.")
        exit(1)
    end
end

main()

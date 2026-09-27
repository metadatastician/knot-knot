# SPDX-License-Identifier: MPL-2.0
# conway-suite.jl — minimal working prototype for the Conway suite
# (src/tangle.jl + src/skein.jl).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/conway-suite.jl`; it exits
# non-zero unless every pinned assertion holds.
#
# What this proves: the rational-tangle calculus is exact rational arithmetic
# (horizontal sum adds, vertical product multiplies, rotation is -1/f, mirror
# is -f), Conway notation round-trips through the Euclidean algorithm, the
# numerator closure N(p/q) is a knot exactly when p is odd with determinant
# |p|, and the Conway skein relation ties the invariants together:
#
#     grad(L+) - grad(L-) = z * grad(L0)
#
# It also shows the cheap determinant bridge: for a knot the Conway
# polynomial has only even powers, so evaluating it at z = 2i (the point where
# z^2 = t + t^-1 - 2 hits t = -1) gives grad(2i) = sum_i c_i (-4)^i, a real
# number whose absolute value is the determinant.

const FAILURES = Ref(0)

function check(label, got, want)
    if got == want
        println("  PASS  $label")
    else
        FAILURES[] += 1
        println("  FAIL  $label: got $got, want $want")
    end
end

# --- rational tangles -------------------------------------------------------

"""A rational tangle, classified by its fraction in Q (0 meaning the zero tangle)."""
struct RationalTangle
    fraction::Rational{Int}
end

"""Horizontal sum: `frac(T + S) = frac(T) + frac(S)`."""
Base.:+(t::RationalTangle, s::RationalTangle) = RationalTangle(t.fraction + s.fraction)

"""Vertical product: `frac(T * S) = frac(T) * frac(S)`."""
Base.:*(t::RationalTangle, s::RationalTangle) = RationalTangle(t.fraction * s.fraction)

"""Mirror tangle: negates the fraction."""
mirror(t::RationalTangle) = RationalTangle(-t.fraction)

"""Quarter-turn rotation: sends the fraction `f` to `-1/f`."""
function rotate(t::RationalTangle)
    t.fraction == 0 && error("rotation of the zero tangle has infinite fraction")
    RationalTangle(-1 // t.fraction)
end

"""
    continued_fraction(terms) -> Rational{Int}

Evaluate `[a1, ..., an]` right-to-left:
`a1 + 1/(a2 + 1/(... + 1/an))`. The empty list is the zero tangle.
"""
function continued_fraction(terms::Vector{Int})
    isempty(terms) && return 0 // 1
    acc = Rational{Int}(terms[end])
    for k in (length(terms) - 1):-1:1
        acc == 0 && error("continued fraction has a zero tail value")
        acc = terms[k] + 1 // acc
    end
    acc
end

"""Recover a Conway notation list for `p//q` by the Euclidean algorithm."""
function conway_notation(fraction::Rational{Int})
    p, q = Int(numerator(fraction)), Int(denominator(fraction))
    q == 0 && error("infinite fraction has no finite Conway notation")
    terms = Int[]
    while q != 0
        a = fld(p, q)
        push!(terms, a)
        p, q = q, p - a * q
    end
    terms
end

RationalTangle(terms::Vector{Int}) = RationalTangle(continued_fraction(terms))

"""Determinant of the 2-bridge link `N(p/q)`: the absolute numerator `|p|`."""
two_bridge_det(fraction::Rational{Int}) = abs(Int(numerator(fraction)))

"""True when `N(p/q)` is a knot: `p` odd. Even numerators give 2-component links."""
is_knot_fraction(fraction::Rational{Int}) = isodd(numerator(fraction))

# --- exact Conway polynomials ----------------------------------------------

# A Conway polynomial in z, as exponent -> Rational{Int} coefficient. Zero
# exponents are dropped so equality is structural.

lp_norm(p::Dict{Int,Rational{Int}}) =
    Dict(k => v for (k, v) in p if !iszero(v))

lp_zero() = Dict{Int,Rational{Int}}()
lp_one() = lp_norm(Dict(0 => 1 // 1))
lp_z() = lp_norm(Dict(1 => 1 // 1))

function lp_add(a::Dict{Int,Rational{Int}}, b::Dict{Int,Rational{Int}})
    out = copy(a)
    for (k, v) in b
        out[k] = get(out, k, 0 // 1) + v
    end
    lp_norm(out)
end

function lp_mul(a::Dict{Int,Rational{Int}}, b::Dict{Int,Rational{Int}})
    out = lp_zero()
    for (i, x) in a, (j, y) in b
        out[i + j] = get(out, i + j, 0 // 1) + x * y
    end
    lp_norm(out)
end

lp_neg(a) = lp_norm(Dict(k => -v for (k, v) in a))
lp_sub(a, b) = lp_add(a, lp_neg(b))
lp_scale(a, c::Rational{Int}) = lp_norm(Dict(k => v * c for (k, v) in a))

"""
    det_via_conway(p) -> Int

The determinant from a Conway polynomial. Substituting `z = 2i` is the point
where `z^2 = t + t^-1 - 2` reaches `t = -1`, so `|grad(2i)| = |Delta(-1)|`.
Because a knot's Conway polynomial has only even powers, `(2i)^(2i) = (-4)^i`
is real and the arithmetic stays exact.
"""
function det_via_conway(p::Dict{Int,Rational{Int}})
    any(isodd(k) for k in keys(p)) &&
        error("only knot-type Conway polynomials (even powers only) are supported")
    total = 0 // 1
    for (k, v) in p
        total += v * (-4 // 1)^(k ÷ 2)
    end
    abs(Int(numerator(total)))
end

# --- the pins ---------------------------------------------------------------

function main()
    println("conway-suite prototype")
    println("  (rational tangles, Conway notation, 2-bridge links, the skein relation)")

    println("\n[1] continued fractions, right-to-left")
    check("[1,1,1,1] = 5/3", continued_fraction([1, 1, 1, 1]), 5 // 3)
    check("[2,2] = 5/2", continued_fraction([2, 2]), 5 // 2)
    check("[3] = 3", continued_fraction([3]), 3 // 1)
    check("[] = 0 (the zero tangle)", continued_fraction(Int[]), 0 // 1)
    check("[1,1,2] = 5/3", continued_fraction([1, 1, 2]), 5 // 3)
    check("[1,1,1,1,1] = 8/5", continued_fraction([1, 1, 1, 1, 1]), 8 // 5)

    println("\n[2] Conway notation round-trips")
    for f in (5 // 3, 5 // 2, 3 // 1, 9 // 2, 2 // 1, 7 // 3, 8 // 5)
        check("conway_notation($f) re-evaluates to $f",
              continued_fraction(conway_notation(f)), f)
    end

    println("\n[3] 2-bridge classification")
    check("det N(3/1)", two_bridge_det(3 // 1), 3)
    check("det N(5/3)", two_bridge_det(5 // 3), 5)
    check("det N(5/2)", two_bridge_det(5 // 2), 5)
    check("det N(9/2)", two_bridge_det(9 // 2), 9)
    check("N(3/1) is a knot", is_knot_fraction(3 // 1), true)
    check("N(5/3) is a knot", is_knot_fraction(5 // 3), true)
    check("N(5/2) is a knot", is_knot_fraction(5 // 2), true)
    check("N(9/2) is a knot", is_knot_fraction(9 // 2), true)
    check("N(2/1) is the Hopf link, not a knot", is_knot_fraction(2 // 1), false)

    println("\n[4] tangle calculus identities")
    for f in (5 // 3, 5 // 2, 3 // 1, 2 // 1)
        t = RationalTangle(f)
        check("rotation squared is the identity on $f",
              rotate(rotate(t)).fraction, f)
        check("mirror and rotation commute on $f",
              mirror(rotate(t)).fraction, rotate(mirror(t)).fraction)
        check("rotation four times is the identity on $f",
              rotate(rotate(rotate(rotate(t)))).fraction, f)
        # horizontal sum and vertical product really do add and multiply
        s = RationalTangle(2 // 1)
        check("horizontal sum adds on $f", (t + s).fraction, f + 2 // 1)
        check("vertical product multiplies on $f", (t * s).fraction, f * 2 // 1)
        # the two sums together recover the fraction of the numerator closure
        check("sum of a tangle and its mirror is 0", (t + mirror(t)).fraction, 0 // 1)
    end
    check("Conway notation builds the same tangle as the fraction",
          RationalTangle([1, 1, 1, 1]).fraction, 5 // 3)

    println("\n[5] the Conway skein relation")
    # grad(3_1) - grad(unknot) = z * grad(Hopf), the defining relation.
    grad_trefoil = lp_add(lp_one(), lp_mul(lp_z(), lp_z()))       # 1 + z^2
    grad_unknot = lp_one()
    grad_hopf = lp_z()                                            # z
    check("grad(3_1) - grad(unknot) = z * grad(Hopf)",
          lp_sub(grad_trefoil, grad_unknot), lp_mul(lp_z(), grad_hopf))
    # and the same relation one level up: grad(4_1) - grad(3_1) = -2z^2 = z * (-2z)
    grad_fig8 = lp_add(lp_one(), lp_neg(lp_mul(lp_z(), lp_z())))  # 1 - z^2
    check("grad(4_1) - grad(3_1) = z * (-2z)",
          lp_sub(grad_fig8, grad_trefoil),
          lp_mul(lp_z(), lp_scale(lp_z(), -2 // 1)))

    println("\n[6] determinant bridge |grad(2i)| = |Delta(-1)|")
    check("det(unknot) via Conway", det_via_conway(grad_unknot), 1)
    check("det(trefoil) via Conway", det_via_conway(grad_trefoil), 3)
    check("det(figure-eight) via Conway", det_via_conway(grad_fig8), 5)
    grad_cinquefoil =
        lp_add(lp_one(), lp_add(lp_scale(lp_mul(lp_z(), lp_z()), 3 // 1),
                                lp_mul(lp_mul(lp_z(), lp_z()), lp_mul(lp_z(), lp_z()))))
    check("grad(5_1) = 1 + 3z^2 + z^4",
          grad_cinquefoil, lp_norm(Dict(0 => 1 // 1, 2 => 3 // 1, 4 => 1 // 1)))
    check("det(cinquefoil) via Conway", det_via_conway(grad_cinquefoil), 5)

    println("\n[7] the determinant agrees with the 2-bridge numerator")
    # N(3/1) is the trefoil, N(5/3) and N(5/2) the figure-eight (they are the
    # same knot up to the q -> q^-1 symmetry), N(5/1) the cinquefoil, so the
    # tangle route and the polynomial route must agree.
    check("N(3/1) is the trefoil: determinant matches grad(3_1)",
          two_bridge_det(3 // 1), det_via_conway(grad_trefoil))
    check("N(5/3) is the figure-eight: determinant matches grad(4_1)",
          two_bridge_det(5 // 3), det_via_conway(grad_fig8))
    check("N(5/2) is also the figure-eight: determinant matches grad(4_1)",
          two_bridge_det(5 // 2), det_via_conway(grad_fig8))
    check("N(5/1) is the cinquefoil: determinant matches grad(5_1)",
          two_bridge_det(5 // 1), det_via_conway(grad_cinquefoil))
    # 5/1, 5/2 and 5/3 all have determinant 5 but are not all the same knot,
    # which is exactly why the determinant is only an obstruction.
    check("determinant 5 covers both 4_1 and 5_1",
          (two_bridge_det(5 // 3), two_bridge_det(5 // 1)), (5, 5))

    println()
    if FAILURES[] == 0
        println("ALL PASS — the Conway-suite prototype is sound.")
    else
        println("$(FAILURES[]) FAILURE(S) — the Conway-suite prototype is wrong.")
        exit(1)
    end
end

main()

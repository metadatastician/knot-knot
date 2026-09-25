# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Rational tangles and Conway calculus.
#
# A rational tangle is classified by its fraction p/q in Q ∪ {∞}; the
# calculus is exact rational arithmetic:
#   * horizontal sum   frac(T + S) = frac(T) + frac(S)
#   * vertical product frac(T * S) = frac(T) * frac(S)
#   * rotation         frac(rot T) = -1 / frac(T)
#   * mirror           frac(mir T) = -frac(T)
# The numerator closure of the tangle p/q is the 2-bridge link N(p/q):
# a knot exactly when p is odd, with determinant |p|. Conway notation
# [a1, ..., an] evaluates right-to-left: a1 + 1/(a2 + 1/(... + 1/an)).

"""
    RationalTangle(fraction::Rational{Int})

A rational tangle, classified by its Conway fraction. The fraction `1//0`
is represented by `typemin`-guarded arithmetic nowhere in this package:
infinite fractions are not constructible; use the ∞ tangle's role through
rotation of the zero tangle only when `q != 0` results.
"""
struct RationalTangle
    fraction::Rational{Int}
end

"""
    tangle_fraction(x) -> Rational{Int}

The Conway fraction of a rational tangle.
"""
tangle_fraction(t::RationalTangle)::Rational{Int} = t.fraction

Base.:(+)(a::RationalTangle, b::RationalTangle) = RationalTangle(a.fraction + b.fraction)
Base.:(*)(a::RationalTangle, b::RationalTangle) = RationalTangle(a.fraction * b.fraction)
Base.:(-)(a::RationalTangle) = RationalTangle(-a.fraction)

"""
    mirror(t::RationalTangle)

Mirror tangle: negates the fraction.
"""
mirror(t::RationalTangle) = -t

"""
    rotate(t::RationalTangle)

Quarter-turn rotation: sends the fraction `f` to `-1/f` (Conway's
reciprocal-with-sign rule).
"""
function rotate(t::RationalTangle)::RationalTangle
    t.fraction == 0 && error("rotation of the zero tangle has infinite fraction")
    RationalTangle(-1 // t.fraction)
end

"""
    continued_fraction(terms::Vector{Int}) -> Rational{Int}

Evaluate the continued fraction `[a1, ..., an]` right-to-left:
``a_1 + \\cfrac{1}{a_2 + \\cfrac{1}{\\ddots + 1/a_n}}``.
"""
function continued_fraction(terms::Vector{Int})::Rational{Int}
    isempty(terms) && return 0 // 1
    acc = Rational{Int}(terms[end])
    for k in (length(terms) - 1):-1:1
        acc == 0 && error("continued fraction has a zero tail value")
        acc = terms[k] + 1 // acc
    end
    acc
end

"""
    conway_notation(fraction::Rational{Int}) -> Vector{Int}

Recover a Conway continued-fraction list `[a1, ..., an]` for `p//q` via the
Euclidean algorithm (floor quotients, right-to-left evaluation convention).
`conway_notation(continued_fraction(x)) == x` round-trips up to the usual
continued-fraction ambiguity.
"""
function conway_notation(fraction::Rational{Int})::Vector{Int}
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

"""
    RationalTangle(terms::Vector{Int})

The rational tangle with Conway notation `[a1, ..., an]`.
"""
RationalTangle(terms::Vector{Int}) = RationalTangle(continued_fraction(terms))

"""
    two_bridge_det(fraction::Rational{Int}) -> Int

Determinant of the 2-bridge link ``N(p/q)``: the absolute numerator ``|p|``.
"""
two_bridge_det(fraction::Rational{Int})::Int = abs(Int(numerator(fraction)))

"""
    is_knot_fraction(fraction::Rational{Int}) -> Bool

True when the numerator closure ``N(p/q)`` is a knot (``p`` odd); even
numerators give two-component links.
"""
is_knot_fraction(fraction::Rational{Int})::Bool = isodd(numerator(fraction))

function Base.show(io::IO, t::RationalTangle)
    f = t.fraction
    print(io, "RationalTangle(", numerator(f), "/", denominator(f), ")")
end

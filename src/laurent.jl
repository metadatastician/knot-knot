# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Laurent polynomials with integer coefficients, keyed by exponent.
# The invariant kernels (Jones, Alexander-Conway) only ever produce
# integer coefficients, so a Dict{Int,Int} core is exact and
# dependency-free. Wrapped in a named struct so all arithmetic is defined
# on KnotKnot's own type (no Base-type overloading).

"""
    LaurentPoly

A univariate Laurent polynomial ``\\sum_e c_e t^e`` with integer
coefficients. Zero coefficients are trimmed eagerly by every arithmetic
operation; equality and display use the trimmed form.
"""
struct LaurentPoly
    terms::Dict{Int,Int}
    LaurentPoly(terms::Dict{Int,Int}) = new(terms)
end

LaurentPoly() = LaurentPoly(Dict{Int,Int}())

function LaurentPoly(pairs::Pair{Int,Int}...)
    d = Dict{Int,Int}()
    for (e, c) in pairs
        d[e] = get(d, e, 0) + c
    end
    LaurentPoly(d)
end

function _lp_trim!(p::LaurentPoly)::LaurentPoly
    filter!(kv -> kv[2] != 0, p.terms)
    p
end

_lp_trimmed(p::LaurentPoly)::LaurentPoly = _lp_trim!(LaurentPoly(copy(p.terms)))

Base.isempty(p::LaurentPoly)::Bool = isempty(p.terms)
Base.keys(p::LaurentPoly) = keys(p.terms)
Base.values(p::LaurentPoly) = values(p.terms)
Base.haskey(p::LaurentPoly, e::Int)::Bool = haskey(p.terms, e)
Base.getindex(p::LaurentPoly, e::Int)::Int = p.terms[e]
Base.get(p::LaurentPoly, e::Int, default::Int)::Int = get(p.terms, e, default)
Base.setindex!(p::LaurentPoly, c::Int, e::Int) = (p.terms[e] = c; p)
Base.iterate(p::LaurentPoly) = iterate(p.terms)
Base.iterate(p::LaurentPoly, st) = iterate(p.terms, st)
Base.iszero(p::LaurentPoly)::Bool = isempty(_lp_trimmed(p))

function Base.:(==)(a::LaurentPoly, b::LaurentPoly)::Bool
    _lp_trimmed(a).terms == _lp_trimmed(b).terms
end

function Base.:(+)(a::LaurentPoly, b::LaurentPoly)::LaurentPoly
    r = copy(a.terms)
    for (e, c) in b
        r[e] = get(r, e, 0) + c
    end
    _lp_trim!(LaurentPoly(r))
end

function Base.:(-)(a::LaurentPoly, b::LaurentPoly)::LaurentPoly
    r = copy(a.terms)
    for (e, c) in b
        r[e] = get(r, e, 0) - c
    end
    _lp_trim!(LaurentPoly(r))
end

Base.:(-)(a::LaurentPoly)::LaurentPoly = LaurentPoly(Dict(e => -c for (e, c) in a))

function Base.:(*)(a::LaurentPoly, b::LaurentPoly)::LaurentPoly
    r = Dict{Int,Int}()
    for (ea, ca) in a, (eb, cb) in b
        r[ea + eb] = get(r, ea + eb, 0) + ca * cb
    end
    _lp_trim!(LaurentPoly(r))
end

Base.:(*)(k::Integer, p::LaurentPoly)::LaurentPoly =
    LaurentPoly(Dict(e => k * c for (e, c) in p))
Base.:(*)(p::LaurentPoly, k::Integer)::LaurentPoly = k * p

"""
    lp_pow(p::LaurentPoly, n::Integer) -> LaurentPoly

Non-negative integer power of a Laurent polynomial. `lp_pow(p, 0)` is the
constant polynomial 1.
"""
function lp_pow(p::LaurentPoly, n::Integer)::LaurentPoly
    n >= 0 || error("lp_pow: negative powers need division, not representable here")
    n == 0 && return LaurentPoly(0 => 1)
    r = LaurentPoly(0 => 1)
    for _ in 1:n
        r = r * p
    end
    r
end

"""
    lp_shift(p::LaurentPoly, k::Integer) -> LaurentPoly

Multiply `p` by ``t^k`` (shift every exponent by `k`).
"""
lp_shift(p::LaurentPoly, k::Integer)::LaurentPoly =
    LaurentPoly(Dict(e + k => c for (e, c) in p))

"""
    lp_trim(p::LaurentPoly) -> LaurentPoly

A copy of `p` with all zero coefficients removed.
"""
lp_trim(p::LaurentPoly)::LaurentPoly = _lp_trimmed(p)

"""
    lp_eval(p::LaurentPoly, x::Number) -> Number

Evaluate `p` at `x`. Negative exponents require a non-zero `x`.
"""
function lp_eval(p::LaurentPoly, x::Number)
    acc = 0 * x
    for (e, c) in p
        acc += c * x^e
    end
    acc
end

"""
    lp_show_string(p::LaurentPoly; var = "t") -> String

Render in conventional mathematical notation, highest exponent first,
e.g. `-t^-4 + t^-3 + t^-1`. The zero polynomial renders as `0`.
"""
function lp_show_string(p::LaurentPoly; var::String = "t")::String
    p = _lp_trimmed(p)
    isempty(p) && return "0"
    terms = String[]
    for e in sort(collect(keys(p)), rev = true)
        c = p[e]
        c == 0 && continue
        mag = abs(c)
        body = if e == 0
            "$(mag)"
        elseif mag == 1
            e == 1 ? var : "$(var)^$(e)"
        else
            e == 1 ? "$(mag)*$(var)" : "$(mag)*$(var)^$(e)"
        end
        push!(terms, c < 0 ? "-$(body)" : body)
    end
    s = join(terms, " + ")
    replace(s, "+ -" => "- ")
end

function Base.show(io::IO, p::LaurentPoly)
    print(io, "LaurentPoly(", lp_show_string(p), ")")
end

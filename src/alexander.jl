# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Alexander-Conway invariants via Fox calculus on the Wirtinger presentation.
#
# For a crossing X[a,b,c,d] (slot a the entering under-arc, slot c the
# exiting one, the over-strand occupying slots b/d) the Fox row is
#   gen(a): +1,   gen(c): -t^k,   gen(over): t^k - 1
# where the exponent k is fixed by the crossing's GEOMETRY, not by its
# sign: k = +1 when the over-strand enters the crossing through slot d,
# and k = -1 when it enters through slot b. The entering slot is read off
# the diagram's slot orientation (a 2-colouring of the slot constraint
# graph). Reversing every strand orientation flips all k simultaneously,
# which changes every row by a unit factor only, so Delta is unaffected by
# the colouring branch. Using the crossing sign instead of this left/right
# datum silently corrupts diagrams whose crossings mix both over-strand
# directions (e.g. the KnotAtlas PD of 6_3); the rule here reproduces the
# published Alexander polynomials of every knot 3_1 .. 7_7.
#
# Any (n-1) x (n-1) minor of the resulting n x n matrix has determinant
# equal to Delta(t) up to units; the result is canonicalised to a symmetric
# window with Delta(1) = +1 (Conway normalisation).

# --- dense polynomial vectors (index i holds the coefficient of t^(i-1)) ---

function _vec_trim(p::Vector{Int})::Vector{Int}
    last = length(p)
    while last > 0 && p[last] == 0
        last -= 1
    end
    p[1:last]
end

function _vec_add(a::Vector{Int}, b::Vector{Int})::Vector{Int}
    r = zeros(Int, max(length(a), length(b)))
    for i in eachindex(a)
        r[i] += a[i]
    end
    for i in eachindex(b)
        r[i] += b[i]
    end
    _vec_trim(r)
end

function _vec_sub(a::Vector{Int}, b::Vector{Int})::Vector{Int}
    r = zeros(Int, max(length(a), length(b)))
    for i in eachindex(a)
        r[i] += a[i]
    end
    for i in eachindex(b)
        r[i] -= b[i]
    end
    _vec_trim(r)
end

function _vec_mul(a::Vector{Int}, b::Vector{Int})::Vector{Int}
    (isempty(a) || isempty(b)) && return Int[]
    r = zeros(Int, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        r[i + j - 1] += a[i] * b[j]
    end
    _vec_trim(r)
end

# Exact division by a previous Bareiss pivot (guaranteed divisible).
function _vec_divexact(a::Vector{Int}, b::Vector{Int})::Vector{Int}
    isempty(b) && error("division by the zero polynomial")
    isempty(a) && return Int[]
    deg_b = length(b) - 1
    q = zeros(Int, length(a) - length(b) + 1)
    work = copy(a)
    lead = b[end]
    for k in (length(q) - 1):-1:0
        isempty(work) && break
        idx = k + deg_b + 1  # 1-based position of the t^(k + deg_b) term
        idx > length(work) && continue  # no term to cancel at this degree
        coeff = div(work[idx], lead)
        coeff == 0 && continue
        q[k + 1] = coeff
        for j in eachindex(b)
            work[k + j] -= coeff * b[j]
        end
        work = _vec_trim(work)
    end
    _vec_trim(q)
end

function _vec_eval(p::Vector{Int}, x::Int)::Int
    acc = 0
    for c in Iterators.reverse(p)
        acc = acc * x + c
    end
    acc
end

# Bareiss fraction-free elimination over Z[t].
function _poly_det_bareiss(M::Matrix{Vector{Int}}, n::Int)::Vector{Int}
    A = Matrix{Vector{Int}}(undef, n, n)
    for i in 1:n, j in 1:n
        A[i, j] = _vec_trim(M[i, j])
    end
    sgn = 1
    prev = [1]
    for k in 1:(n - 1)
        piv = 0
        for r in k:n
            if !isempty(A[r, k])
                piv = r
                break
            end
        end
        piv == 0 && return Int[]
        if piv != k
            for c in 1:n
                A[k, c], A[piv, c] = A[piv, c], A[k, c]
            end
            sgn = -sgn
        end
        for i in (k + 1):n
            for j in (k + 1):n
                num = _vec_sub(_vec_mul(A[k, k], A[i, j]), _vec_mul(A[i, k], A[k, j]))
                A[i, j] = _vec_divexact(num, prev)
            end
            A[i, k] = Int[]
        end
        prev = A[k, k]
    end
    d = A[n, n]
    sgn == -1 ? Int[-c for c in d] : d
end

"""
    alexander_polynomial(d::PlanarDiagram) -> LaurentPoly

The Alexander polynomial ``\\Delta(t)`` of the diagram, computed by Fox
calculus on the Wirtinger presentation and normalised so that the exponent
window is symmetric and ``\\Delta(1) = +1`` (Conway normalisation).

Pinned values: trefoil ``t^{-1} - 1 + t``, figure-eight ``-t^{-1} + 3 - t``.
"""
function alexander_polynomial(d::PlanarDiagram)::LaurentPoly
    n = length(d.crossings)
    n == 0 && return LaurentPoly(0 => 1)
    gen_map, n_gens = wirt_generators(d)
    n_gens == 0 && return LaurentPoly(0 => 1)

    o = slot_orientation(d)
    M = [LaurentPoly() for _ in 1:n, _ in 1:n_gens]
    for (i, c) in enumerate(d.crossings)
        a, _, cc, dd = c.arcs
        ga = gen_map[a]
        gc = gen_map[cc]
        gover = gen_map[dd]
        # k = +1 when the over-strand enters through slot d, -1 through b.
        kap = o[4 * (i - 1) + 4] == 1 ? 1 : -1
        M[i, ga][0] = get(M[i, ga], 0, 0) + 1
        M[i, gc][kap] = get(M[i, gc], kap, 0) - 1
        M[i, gover][kap] = get(M[i, gover], kap, 0) + 1
        M[i, gover][0] = get(M[i, gover], 0, 0) - 1
    end

    minor_size = min(n, n_gens) - 1
    minor_size <= 0 && return LaurentPoly(0 => 1)

    # For a knot every cofactor of the Alexander matrix is +-t^k * Delta,
    # but for a given diagram a particular row/column choice can still be a
    # degenerate representative; try a small set of minors before declaring
    # the determinant (and hence the polynomial) zero.
    row_choices = unique([1:minor_size, (n - minor_size + 1):n])
    col_choices = unique([1:minor_size, (n_gens - minor_size + 1):n_gens])
    det_coeffs = Int[]
    shift = 0
    for rows in row_choices, cols in col_choices
        min_exp = 0
        max_exp = 0
        for i in rows, j in cols
            for e in keys(M[i, j])
                min_exp = min(min_exp, e)
                max_exp = max(max_exp, e)
            end
        end
        s = -min_exp
        poly_size = max_exp - min_exp + 1
        PM = Matrix{Vector{Int}}(undef, minor_size, minor_size)
        for (ri, i) in enumerate(rows), (cj, j) in enumerate(cols)
            coeffs = zeros(Int, poly_size)
            for (e, c) in M[i, j]
                c != 0 && (coeffs[e + s + 1] += c)
            end
            PM[ri, cj] = coeffs
        end
        cand = _poly_det_bareiss(PM, minor_size)
        if !isempty(cand)
            det_coeffs = cand
            shift = s
            break
        end
    end
    isempty(det_coeffs) && return LaurentPoly()  # zero polynomial

    overall_shift = shift * minor_size
    poly = LaurentPoly()
    for (k, c) in enumerate(det_coeffs)
        c != 0 && (poly[k - 1 - overall_shift] = c)
    end
    isempty(poly) && return LaurentPoly(0 => 1)

    # Canonical window: symmetric about 0 (even span) or 1/2 (odd span).
    min_e = minimum(keys(poly))
    max_e = maximum(keys(poly))
    span = max_e - min_e
    target_sum = isodd(span) ? 1 : 0
    center_shift = div(target_sum - (min_e + max_e), 2)
    if center_shift != 0
        poly = lp_shift(poly, center_shift)
    end

    # Conway normalisation: Delta(1) = +1 for knots; fall back to a positive
    # leading coefficient for links where Delta(1) = 0.
    delta_at_1 = sum(values(poly))
    flip = if delta_at_1 != 0
        delta_at_1 < 0
    else
        poly[maximum(keys(poly))] < 0
    end
    if flip
        for e in keys(poly)
            poly[e] = -poly[e]
        end
    end
    lp_trim(poly)
end

"""
    determinant(d::PlanarDiagram) -> Int

The knot determinant ``|\\Delta(-1)|``: a positive odd integer for knots.
"""
function determinant(d::PlanarDiagram)::Int
    alex = alexander_polynomial(d)
    val = sum(c * (-1)^e for (e, c) in alex)
    abs(val) == 0 ? 1 : abs(val)
end

"""
    conway_polynomial(d::PlanarDiagram) -> LaurentPoly

The Conway polynomial ``\\nabla(z)`` obtained from ``\\Delta`` by the
substitution ``\\Delta(t) = \\nabla(t^{1/2} - t^{-1/2})``. Even-component
links (odd-span ``\\Delta``) take the antisymmetric branch so that, for
example, the positive Hopf link gives ``\\nabla = z``.
"""
function conway_polynomial(d::PlanarDiagram)::LaurentPoly
    alex = alexander_polynomial(d)
    isempty(alex) && return LaurentPoly()  # zero polynomial (e.g. split links)

    if isodd(maximum(keys(alex)) - minimum(keys(alex)))
        result = _conway_from_antisymmetric(alex)
        # Alexander data fixes an even-component link's Conway potential only
        # up to overall sign; pin the sign with the linking number, whose
        # value is the coefficient of z.
        if !isempty(result) && length(traverse_components(d)) == 2
            lk = linking_number(d)
            c1 = get(result, 1, 0)
            if lk != 0 && c1 != 0 && sign(lk) != sign(c1)
                result = -result
            end
        end
        return result
    end

    max_exp = maximum(abs(e) for e in keys(alex))
    a0 = get(alex, 0, 0)

    # Chebyshev recurrence for S_k(w) = t^k + t^-k with w = z^2 + 2.
    S = Vector{Vector{Int}}(undef, max_exp + 1)
    S[1] = [2]
    if max_exp >= 1
        S[2] = [0, 1]
    end
    for k in 2:max_exp
        wp = vcat([0], S[k])
        r = zeros(Int, max(length(wp), length(S[k - 1])))
        for i in eachindex(wp)
            r[i] += wp[i]
        end
        for i in eachindex(S[k - 1])
            r[i] -= S[k - 1][i]
        end
        S[k + 1] = r
    end

    combined_w = [a0]
    for k in 1:max_exp
        ak = get(alex, k, 0)
        ak == 0 && continue
        combined_w = _vec_add(combined_w, [ak * c for c in S[k + 1]])
    end

    Pz = Vector{Vector{Int}}(undef, length(combined_w))
    Pz[1] = [1]
    for k in 2:length(combined_w)
        Pz[k] = _vec_mul(Pz[k - 1], [2, 1])
    end

    nabla = Int[]
    for (k, ck) in enumerate(combined_w)
        ck == 0 && continue
        nabla = _vec_add(nabla, [ck * c for c in Pz[k]])
    end

    result = LaurentPoly()
    for (k, c) in enumerate(nabla)
        c != 0 && (result[2 * (k - 1)] = c)
    end
    isempty(result) ? LaurentPoly(0 => 1) : result
end

function _conway_from_antisymmetric(alex::LaurentPoly)::LaurentPoly
    max_e = maximum(keys(alex))
    m_max = 2 * max_e - 1
    m_max < 1 && return LaurentPoly(0 => 0)
    w = [2, 0, 1]
    D = Dict{Int,Vector{Int}}()
    D[1] = [0, 1]
    if m_max >= 3
        D[3] = _vec_mul(D[1], [3, 0, 1])
    end
    for m in 5:2:m_max
        D[m] = _vec_sub(_vec_mul(w, D[m - 2]), D[m - 4])
    end
    nabla = Int[]
    for m in 1:2:m_max
        c = get(alex, div(m + 1, 2), 0)
        c == 0 && continue
        nabla = _vec_add(nabla, [c * x for x in D[m]])
    end
    result = LaurentPoly()
    for (i, c) in enumerate(nabla)
        c != 0 && (result[i - 1] = c)
    end
    isempty(result) ? LaurentPoly(0 => 0) : result
end

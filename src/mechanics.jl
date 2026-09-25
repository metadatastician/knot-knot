# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Applied knot mechanics: elastic rods, capstan friction, bend-stress
# superposition, and rope-efficiency prediction.
#
# Unit convention (engineering): lengths in mm, force in N, modulus E in
# MPa (1 GPa = 1000 MPa), stress in MPa (= N/mm^2). The theory-practice
# bridge formula superposing axial and bending stress at the strand surface:
#
#     sigma_total = F_local / A  +  E * kappa * d / 2
#
# where A = pi (d/2)^2 and kappa = 1/R is the local curvature. Inside the
# knot core the carried tension decays by capstan friction,
# F(s) = F0 * exp(-mu * theta(s)), with theta the accumulated contact angle.

"""
    Rod(points, diameter, youngs_mpa, uts_mpa, friction)

A discretised elastic rod centreline.

# Fields
- `points::Vector{NTuple{3,Float64}}`: centreline samples in mm;
- `diameter::Float64`: strand diameter in mm;
- `youngs_mpa::Float64`: Young's modulus in MPa;
- `uts_mpa::Float64`: ultimate tensile strength in MPa;
- `friction::Float64`: strand-on-strand friction coefficient mu.
"""
struct Rod
    points::Vector{NTuple{3,Float64}}
    diameter::Float64
    youngs_mpa::Float64
    uts_mpa::Float64
    friction::Float64
end

function Rod(
    points::Vector{NTuple{3,Float64}};
    diameter::Float64 = 2.0,
    youngs_mpa::Float64 = 3500.0,
    uts_mpa::Float64 = 850.0,
    friction::Float64 = 0.35,
)::Rod
    length(points) >= 3 || error("a rod needs at least three centreline samples")
    Rod(points, diameter, youngs_mpa, uts_mpa, friction)
end

_sub(u::NTuple{3,Float64}, v::NTuple{3,Float64}) = (u[1] - v[1], u[2] - v[2], u[3] - v[3])
_add(u::NTuple{3,Float64}, v::NTuple{3,Float64}) = (u[1] + v[1], u[2] + v[2], u[3] + v[3])
_scale(k::Float64, u::NTuple{3,Float64}) = (k * u[1], k * u[2], k * u[3])
_dot(u::NTuple{3,Float64}, v::NTuple{3,Float64}) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3]
_cross(u::NTuple{3,Float64}, v::NTuple{3,Float64}) =
    (u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1])
_norm(u::NTuple{3,Float64}) = sqrt(_dot(u, u))

"""
    curvature_profile(rod::Rod) -> Vector{Float64}

Discrete curvature ``\\kappa`` (1/mm) at every sample by centred finite
differences: ``\\kappa = |v \\times a| / |v|^3``. End samples mirror their
neighbours (they are usually straight tails).
"""
function curvature_profile(rod::Rod)::Vector{Float64}
    pts = rod.points
    n = length(pts)
    kappa = zeros(Float64, n)
    for i in 2:(n - 1)
        v = _scale(0.5, _sub(pts[i + 1], pts[i - 1]))
        a = _sub(_add(pts[i + 1], pts[i - 1]), _scale(2.0, pts[i]))
        vlen = _norm(v)
        kappa[i] = vlen > 1e-12 ? _norm(_cross(v, a)) / vlen^3 : 0.0
    end
    kappa[1] = kappa[2]
    kappa[n] = kappa[n - 1]
    kappa
end

"""
    segment_lengths(rod::Rod) -> Vector{Float64}

Length of each of the `n-1` centreline segments.
"""
function segment_lengths(rod::Rod)::Vector{Float64}
    [_norm(_sub(rod.points[i + 1], rod.points[i])) for i in 1:(length(rod.points) - 1)]
end

"""
    total_length(rod::Rod) -> Float64

Polyline length of the centreline in mm.
"""
total_length(rod::Rod)::Float64 = sum(segment_lengths(rod))

"""
    bending_energy(rod::Rod) -> Float64

Euler-Bernoulli bending energy ``\\int \\tfrac{1}{2} E I \\kappa^2 \\, ds``
in N*mm, with ``I = \\pi d^4 / 64`` (solid circular strand).
"""
function bending_energy(rod::Rod)::Float64
    kappa = curvature_profile(rod)
    segs = segment_lengths(rod)
    inertia = pi * rod.diameter^4 / 64.0
    acc = 0.0
    for i in eachindex(segs)
        kbar = 0.5 * (kappa[i] + kappa[i + 1])
        acc += 0.5 * rod.youngs_mpa * inertia * kbar^2 * segs[i]
    end
    acc
end

"""
    capstan_tensions(rod::Rod, applied_force::Float64; threshold = 0.5) -> Vector{Float64}

Tension carried along the rod when `applied_force` N is pulled at the tail.
Samples whose curvature is at least `threshold` times the maximum are the
contact core; the tension decays there by the capstan law with the
accumulated tangent turning angle ``\\theta``: ``F(s) = F_0 e^{-\\mu \\theta(s)}``.
"""
function capstan_tensions(
    rod::Rod,
    applied_force::Float64;
    threshold::Float64 = 0.5,
)::Vector{Float64}
    kappa = curvature_profile(rod)
    kmax = maximum(kappa)
    n = length(rod.points)
    in_core = kmax > 0 ? [k >= threshold * kmax for k in kappa] : falses(n)
    tensions = zeros(Float64, n)
    theta = 0.0
    tensions[1] = applied_force
    for i in 2:n
        if in_core[i] && in_core[i - 1]
            t1 = _sub(rod.points[i], rod.points[i - 1])
            l1 = _norm(t1)
            if i + 1 <= n
                t2 = _sub(rod.points[i + 1], rod.points[i])
                l2 = _norm(t2)
                if l1 > 1e-12 && l2 > 1e-12
                    c = clamp(_dot(t1, t2) / (l1 * l2), -1.0, 1.0)
                    theta += acos(c)
                end
            end
        end
        tensions[i] = applied_force * exp(-rod.friction * theta)
    end
    tensions
end

"""
    stress_profile(rod::Rod, applied_force::Float64) -> NamedTuple

Superposed stress state along the rod in MPa:
`axial` (local capstan tension over the cross-section), `bending`
(``E \\kappa d / 2`` at the outer fibre) and `total`.
"""
function stress_profile(rod::Rod, applied_force::Float64)
    area = pi * (rod.diameter / 2.0)^2
    kappa = curvature_profile(rod)
    tensions = capstan_tensions(rod, applied_force)
    axial = [t / area for t in tensions]
    bending = [rod.youngs_mpa * k * rod.diameter / 2.0 for k in kappa]
    (
        axial = axial,
        bending = bending,
        total = [axial[i] + bending[i] for i in eachindex(axial)],
    )
end

"""
    safety_factor(rod::Rod, applied_force::Float64) -> Float64

UTS divided by the peak total stress. Values below 1 mean the rod fails at
this load.
"""
function safety_factor(rod::Rod, applied_force::Float64)::Float64
    sigma = stress_profile(rod, applied_force)
    peak = maximum(sigma.total)
    peak <= 0 && return Inf
    rod.uts_mpa / peak
end

"""
    breaking_force(rod::Rod) -> Float64

Predicted tensile failure load in N: the largest `F` whose peak superposed
stress reaches UTS. Bending stress is load-independent; the axial part
scales linearly with `F`, so the solve is exact per sample.
"""
function breaking_force(rod::Rod)::Float64
    area = pi * (rod.diameter / 2.0)^2
    kappa = curvature_profile(rod)
    unit_tensions = capstan_tensions(rod, 1.0)
    bending = [rod.youngs_mpa * k * rod.diameter / 2.0 for k in kappa]
    best = Inf
    for i in eachindex(kappa)
        alpha = unit_tensions[i] / area
        alpha <= 0 && continue
        cand = (rod.uts_mpa - bending[i]) / alpha
        cand < best && (best = cand)
    end
    max(0.0, best)
end

"""
    knot_efficiency(rod::Rod) -> Float64

Predicted knot strength efficiency: breaking force of the knotted rod over
the straight-strand strength ``\\sigma_{\\mathrm{uts}} A``. Empirical rope
knots sit near 0.60-0.75; this model reproduces that range when the knot
geometry carries a tight nipping turn.
"""
function knot_efficiency(rod::Rod)::Float64
    area = pi * (rod.diameter / 2.0)^2
    breaking_force(rod) / (rod.uts_mpa * area)
end

"""
    writhe_of(points::Vector{NTuple{3,Float64}}) -> Float64

Discrete writhe of a closed polygon by Gauss-integral quadrature:
``\\mathrm{Wr} = \\frac{1}{4\\pi} \\oint\\oint
\\frac{(\\dot r \\times \\dot r') \\cdot (r - r')}{|r - r'|^3} ds\\, ds'``,
with adjacent-sample pairs omitted (their contribution is the local,
integrable part). Planar curves give 0 by symmetry.
"""
function writhe_of(points::Vector{NTuple{3,Float64}})::Float64
    n = length(points)
    n >= 4 || return 0.0
    ds = total_length(Rod(points)) / n
    acc = 0.0
    for i in 1:n, j in 1:n
        gap = abs(i - j)
        (gap <= 1 || gap >= n - 1) && continue
        ri = points[i]
        rj = points[j]
        ti = _sub(points[mod1(i + 1, n)], points[mod1(i - 1, n)])
        tj = _sub(points[mod1(j + 1, n)], points[mod1(j - 1, n)])
        d = _sub(ri, rj)
        r = _norm(d)
        r <= 1e-12 && continue
        acc += _dot(_cross(ti, tj), d) / r^3
    end
    acc * ds * ds / (4.0 * pi)
end

"""
    linking_of(curve_a, curve_b) -> Float64

Linking number of two closed polygons via the Gauss double integral.
A positive Hopf pair gives +1 (orientation-dependent).
"""
function linking_of(
    curve_a::Vector{NTuple{3,Float64}},
    curve_b::Vector{NTuple{3,Float64}},
)::Float64
    na, nb = length(curve_a), length(curve_b)
    (na >= 3 && nb >= 3) || return 0.0
    dsa = total_length(Rod(curve_a)) / na
    dsb = total_length(Rod(curve_b)) / nb
    acc = 0.0
    for i in 1:na, j in 1:nb
        ri = curve_a[i]
        rj = curve_b[j]
        ti = _sub(curve_a[mod1(i + 1, na)], curve_a[mod1(i - 1, na)])
        tj = _sub(curve_b[mod1(j + 1, nb)], curve_b[mod1(j - 1, nb)])
        d = _sub(ri, rj)
        r = _norm(d)
        r <= 1e-12 && continue
        acc += _dot(_cross(ti, tj), d) / r^3
    end
    acc * dsa * dsb / (4.0 * pi)
end

"""
    calugareanu_check(points; twist_turns, ribbon_width) -> NamedTuple

Numerical Calugareanu-White-Fuller identity check for a closed centreline:
builds a ribbon whose material frame rotates `twist_turns` full turns along
the curve, then returns `(linking, writhe, twist, residual)` where
`residual = Lk - (Tw + Wr)` should be near zero. The identity underpins
DNA topology (linking deficit = writhe) in both protistology and
biomechanics applications.
"""
function calugareanu_check(
    points::Vector{NTuple{3,Float64}};
    twist_turns::Float64 = 1.0,
    ribbon_width::Float64 = 0.05,
)
    n = length(points)
    n >= 4 || error("calugareanu_check needs a closed curve with >= 4 samples")
    # A crude but consistent normal field: radial from the centroid.
    cx = sum(p[1] for p in points) / n
    cy = sum(p[2] for p in points) / n
    cz = sum(p[3] for p in points) / n
    edge = Vector{NTuple{3,Float64}}(undef, n)
    for i in 1:n
        p = points[i]
        u = (p[1] - cx, p[2] - cy, p[3] - cz)
        un = _norm(u)
        un <= 1e-12 && (u = (1.0, 0.0, 0.0); un = 1.0)
        u = _scale(1.0 / un, u)
        phi = 2.0 * pi * twist_turns * (i - 1) / n
        # rotate u about the local tangent by phi via two-frame blend
        t = _sub(points[mod1(i + 1, n)], points[mod1(i - 1, n)])
        tn = _norm(t)
        t = tn > 1e-12 ? _scale(1.0 / tn, t) : (0.0, 0.0, 1.0)
        w = _cross(t, u)
        wn = _norm(w)
        w = wn > 1e-12 ? _scale(1.0 / wn, w) : (0.0, 1.0, 0.0)
        dir = _add(_scale(cos(phi), u), _scale(sin(phi), w))
        edge[i] = _add(p, _scale(ribbon_width, dir))
    end
    (
        linking = linking_of(points, edge),
        writhe = writhe_of(points),
        twist = twist_turns,
        residual = linking_of(points, edge) - (twist_turns + writhe_of(points)),
    )
end

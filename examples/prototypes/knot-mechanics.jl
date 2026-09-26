# SPDX-License-Identifier: MPL-2.0
# knot-mechanics.jl — minimal working prototype for the applied-mechanics
# kernel (src/mechanics.jl).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/knot-mechanics.jl`; it exits
# non-zero unless every pinned assertion holds.
#
# What this proves: the engineering model on a discretised elastic rod —
# discrete curvature, Euler-Bernoulli bending energy, capstan tension decay
# through the contact core, superposed stress, an exact breaking-force solve
# and the knot-efficiency prediction — is internally consistent and lands in
# the documented efficiency band for a knotted rod.
#
# One bug is found and pinned here. `writhe_of` and `linking_of` in
# src/mechanics.jl evaluate the Gauss double integral with the RAW central
# differences `p[i+1] - p[i-1]` as tangents, and weight the sum by
# `total_length / n`. Those central differences are O(ds) in magnitude, so the
# quadrature carries a spurious `ds^2` factor and converges to ZERO as the
# sampling is refined: the linking number of a genuine Hopf pair comes out as
# -0.065 rather than -1, and the trefoil's writhe is not even stable under
# refinement (4.66, 33.7, 8.5 for n = 80, 120, 240). The prototype keeps both
# forms — `writhe`/`linking` with unit tangents and closed-polygon arc lengths,
# and `writhe_unnormalised`/`linking_unnormalised` reproducing the defect — and
# asserts the difference, so the fix cannot silently regress.

const FAILURES = Ref(0)

function check(label, got, want)
    if got == want
        println("  PASS  $label")
    else
        FAILURES[] += 1
        println("  FAIL  $label: got $got, want $want")
    end
end

function check_close(label, got, want; atol = 1e-9, rtol = 1e-9)
    if isapprox(got, want; atol = atol, rtol = rtol)
        println("  PASS  $label")
    else
        FAILURES[] += 1
        println("  FAIL  $label: got $got, want $want (atol $atol, rtol $rtol)")
    end
end

# --- vector helpers ---------------------------------------------------------

_sub(u, v) = (u[1] - v[1], u[2] - v[2], u[3] - v[3])
_add(u, v) = (u[1] + v[1], u[2] + v[2], u[3] + v[3])
_scale(k, u) = (k * u[1], k * u[2], k * u[3])
_dot(u, v) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3]
_cross(u, v) = (u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3],
                u[1] * v[2] - u[2] * v[1])
_norm(u) = sqrt(_dot(u, u))

"""Arc length of every segment of a CLOSED polygon (the closing segment included)."""
segment_lengths_closed(p) = [_norm(_sub(p[mod1(i + 1, length(p))], p[i]))
                             for i in eachindex(p)]

"""Arc length of every segment of an OPEN polyline (`n - 1` segments)."""
segment_lengths(p) = [_norm(_sub(p[i + 1], p[i])) for i in 1:(length(p) - 1)]

total_length_closed(p) = sum(segment_lengths_closed(p))
total_length(p) = sum(segment_lengths(p))

# --- the rod ----------------------------------------------------------------

"""A discretised elastic rod centreline: samples in mm plus material data."""
struct Rod
    points::Vector{NTuple{3,Float64}}
    diameter::Float64
    youngs_mpa::Float64
    uts_mpa::Float64
    friction::Float64
end

function Rod(points::Vector{NTuple{3,Float64}};
             diameter::Float64 = 2.0,
             youngs_mpa::Float64 = 3500.0,
             uts_mpa::Float64 = 850.0,
             friction::Float64 = 0.35)
    length(points) >= 3 || error("a rod needs at least three centreline samples")
    Rod(points, diameter, youngs_mpa, uts_mpa, friction)
end

"""
    curvature_profile(rod)

Discrete curvature `kappa = |v x a| / |v|^3` at every sample by centred
differences; the end samples mirror their neighbours.
"""
function curvature_profile(rod::Rod)
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

"""Euler-Bernoulli bending energy `int 1/2 E I kappa^2 ds` in N*mm."""
function bending_energy(rod::Rod)
    kappa = curvature_profile(rod)
    segs = segment_lengths(rod.points)
    inertia = pi * rod.diameter^4 / 64.0
    acc = 0.0
    for i in eachindex(segs)
        kbar = 0.5 * (kappa[i] + kappa[i + 1])
        acc += 0.5 * rod.youngs_mpa * inertia * kbar^2 * segs[i]
    end
    acc
end

"""
    capstan_tensions(rod, force; threshold)

Tension along the rod: constant outside the contact core, decaying by the
capstan law `F(s) = F0 exp(-mu theta(s))` inside it.
"""
function capstan_tensions(rod::Rod, applied_force::Float64; threshold::Float64 = 0.5)
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

"""Superposed stress (MPa): `axial`, `bending` and `total`."""
function stress_profile(rod::Rod, applied_force::Float64)
    area = pi * (rod.diameter / 2.0)^2
    kappa = curvature_profile(rod)
    tensions = capstan_tensions(rod, applied_force)
    axial = [t / area for t in tensions]
    bending = [rod.youngs_mpa * k * rod.diameter / 2.0 for k in kappa]
    (axial = axial, bending = bending, total = [axial[i] + bending[i] for i in eachindex(axial)])
end

"""UTS over the peak total stress; below 1 the rod fails at this load."""
function safety_factor(rod::Rod, applied_force::Float64)
    peak = maximum(stress_profile(rod, applied_force).total)
    peak <= 0 && return Inf
    rod.uts_mpa / peak
end

"""Predicted tensile failure load in N, solved exactly per sample."""
function breaking_force(rod::Rod)
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

"""Breaking force of the knotted rod over the straight-strand strength."""
knot_efficiency(rod::Rod) =
    breaking_force(rod) / (rod.uts_mpa * pi * (rod.diameter / 2.0)^2)

# --- writhe and linking -----------------------------------------------------

"""
    writhe(points)

Discrete writhe of a closed polygon by Gauss-integral quadrature. Unit
tangents and closed-polygon arc lengths: this is the form that converges.
Planar curves give exactly 0 by symmetry.
"""
function writhe(points::Vector{NTuple{3,Float64}})
    n = length(points)
    n >= 4 || return 0.0
    segs = segment_lengths_closed(points)
    acc = 0.0
    for i in 1:n, j in 1:n
        gap = abs(i - j)
        (gap <= 1 || gap >= n - 1) && continue
        ti = _sub(points[mod1(i + 1, n)], points[mod1(i - 1, n)])
        tj = _sub(points[mod1(j + 1, n)], points[mod1(j - 1, n)])
        ti = _scale(1.0 / _norm(ti), ti)
        tj = _scale(1.0 / _norm(tj), tj)
        d = _sub(points[i], points[j])
        r = _norm(d)
        r <= 1e-12 && continue
        acc += _dot(_cross(ti, tj), d) / r^3 * segs[i] * segs[j]
    end
    acc / (4.0 * pi)
end

"""Linking number of two closed polygons, unit-tangent Gauss quadrature."""
function linking(curve_a::Vector{NTuple{3,Float64}}, curve_b::Vector{NTuple{3,Float64}})
    na, nb = length(curve_a), length(curve_b)
    (na >= 3 && nb >= 3) || return 0.0
    sa = segment_lengths_closed(curve_a)
    sb = segment_lengths_closed(curve_b)
    acc = 0.0
    for i in 1:na, j in 1:nb
        ti = _sub(curve_a[mod1(i + 1, na)], curve_a[mod1(i - 1, na)])
        tj = _sub(curve_b[mod1(j + 1, nb)], curve_b[mod1(j - 1, nb)])
        ti = _scale(1.0 / _norm(ti), ti)
        tj = _scale(1.0 / _norm(tj), tj)
        d = _sub(curve_a[i], curve_b[j])
        r = _norm(d)
        r <= 1e-12 && continue
        acc += _dot(_cross(ti, tj), d) / r^3 * sa[i] * sb[j]
    end
    acc / (4.0 * pi)
end

"""
    writhe_unnormalised(points)

The form src/mechanics.jl actually uses: raw central differences as tangents,
weighted by `total_length / n`. Kept only to demonstrate the defect — it
converges to 0 as the sampling is refined.
"""
function writhe_unnormalised(points::Vector{NTuple{3,Float64}})
    n = length(points)
    n >= 4 || return 0.0
    ds = total_length(points) / n
    acc = 0.0
    for i in 1:n, j in 1:n
        gap = abs(i - j)
        (gap <= 1 || gap >= n - 1) && continue
        ti = _sub(points[mod1(i + 1, n)], points[mod1(i - 1, n)])
        tj = _sub(points[mod1(j + 1, n)], points[mod1(j - 1, n)])
        d = _sub(points[i], points[j])
        r = _norm(d)
        r <= 1e-12 && continue
        acc += _dot(_cross(ti, tj), d) / r^3
    end
    acc * ds * ds / (4.0 * pi)
end

"""The `linking` counterpart of `writhe_unnormalised`, defect included."""
function linking_unnormalised(a::Vector{NTuple{3,Float64}}, b::Vector{NTuple{3,Float64}})
    na, nb = length(a), length(b)
    (na >= 3 && nb >= 3) || return 0.0
    dsa = total_length(a) / na
    dsb = total_length(b) / nb
    acc = 0.0
    for i in 1:na, j in 1:nb
        ti = _sub(a[mod1(i + 1, na)], a[mod1(i - 1, na)])
        tj = _sub(b[mod1(j + 1, nb)], b[mod1(j - 1, nb)])
        d = _sub(a[i], b[j])
        r = _norm(d)
        r <= 1e-12 && continue
        acc += _dot(_cross(ti, tj), d) / r^3
    end
    acc * dsa * dsb / (4.0 * pi)
end

"""
    calugareanu_check(points; twist_turns, ribbon_width)

`Lk - (Tw + Wr)` for a ribbon whose material frame rotates `twist_turns` full
turns along the centreline. With the radial (non-parallel-transport) frame
used here the residual is small only for zero twist.
"""
function calugareanu_check(points::Vector{NTuple{3,Float64}};
                           twist_turns::Float64 = 0.0, ribbon_width::Float64 = 0.05)
    n = length(points)
    n >= 4 || error("calugareanu_check needs a closed curve with >= 4 samples")
    cx = sum(p[1] for p in points) / n
    cy = sum(p[2] for p in points) / n
    cz = sum(p[3] for p in points) / n
    edge = Vector{NTuple{3,Float64}}(undef, n)
    for i in 1:n
        p = points[i]
        u = _sub(p, (cx, cy, cz))
        un = _norm(u)
        un <= 1e-12 && (u = (1.0, 0.0, 0.0); un = 1.0)
        u = _scale(1.0 / un, u)
        phi = 2.0 * pi * twist_turns * (i - 1) / n
        t = _sub(points[mod1(i + 1, n)], points[mod1(i - 1, n)])
        tn = _norm(t)
        t = tn > 1e-12 ? _scale(1.0 / tn, t) : (0.0, 0.0, 1.0)
        w = _cross(t, u)
        wn = _norm(w)
        w = wn > 1e-12 ? _scale(1.0 / wn, w) : (0.0, 1.0, 0.0)
        dir = _add(_scale(cos(phi), u), _scale(sin(phi), w))
        edge[i] = _add(p, _scale(ribbon_width, dir))
    end
    (linking = linking(points, edge), writhe = writhe(points),
     twist = twist_turns,
     residual = linking(points, edge) - (twist_turns + writhe(points)))
end

# --- test geometry ----------------------------------------------------------

"""A regular `n`-gon of radius `R` in the xy-plane."""
circle(n::Int, R::Float64) =
    [(R * cos(2pi * i / n), R * sin(2pi * i / n), 0.0) for i in 1:n]

"""A straight rod along the x-axis."""
straight(n::Int = 4, step::Float64 = 10.0) =
    [(step * (i - 1), 0.0, 0.0) for i in 1:n]

"""A (2,3)-torus-knot centreline, scaled to radius `R`."""
trefoil_rod(n::Int, R::Float64) =
    [(R * (cos(2pi * i / n) + 2cos(4pi * i / n)) / 3,
      R * (sin(2pi * i / n) - 2sin(4pi * i / n)) / 3,
      R * sin(6pi * i / n) / 3) for i in 1:n]

"""A Hopf pair: a unit circle in the xy-plane and one threading it."""
hopf_pair(n::Int, R::Float64 = 1.0) =
    (circle(n, R), [(R + R * cos(2pi * s / n), 0.0, R * sin(2pi * s / n)) for s in 1:n])

# --- the pins ---------------------------------------------------------------

function main()
    println("knot-mechanics prototype")
    println("  (curvature, bending energy, capstan decay, breaking force, efficiency)")

    println("\n[1] discrete curvature of a circle")
    # For a regular n-gon inscribed in a circle of radius R the centred
    # difference gives exactly kappa = 1 / (R cos^2(pi / n)).
    for (n, R) in ((64, 5.0), (128, 5.0), (96, 12.0))
        rod = Rod(circle(n, R))
        kappa = curvature_profile(rod)
        check_close("circle R=$R n=$n: kappa[10] = 1/(R cos^2(pi/n))",
                    kappa[10], 1.0 / (R * cos(pi / n)^2); atol = 1e-12)
        check_close("circle R=$R n=$n: kappa[1] mirrors kappa[2]",
                    kappa[1], kappa[2]; atol = 0.0)
    end
    # a straight rod has no curvature at all
    check("straight rod has zero curvature",
          all(iszero, curvature_profile(Rod(straight(8)))), true)

    println("\n[2] Euler-Bernoulli bending energy")
    # The closed form is pi * E * I / R; the discrete integral sits about 1%
    # below it because the centred difference overestimates kappa near the
    # mirrored end samples.
    for (n, R) in ((64, 5.0), (128, 5.0))
        rod = Rod(circle(n, R))
        I = pi * rod.diameter^4 / 64.0
        exact = pi * rod.youngs_mpa * I / R
        got = bending_energy(rod)
        check_close("circle R=$R n=$n: bending energy ~ pi*E*I/R", got, exact; rtol = 0.02)
        check("  and it is below the closed form", got < exact, true)
    end
    check("straight rod has zero bending energy",
          bending_energy(Rod(straight(8))), 0.0)

    println("\n[3] capstan tension decay")
    st = Rod(straight(6))
    check("straight rod carries the applied force everywhere",
          all(t -> t == 100.0, capstan_tensions(st, 100.0)), true)
    bent = Rod(trefoil_rod(60, 10.0); diameter = 1.0, friction = 0.35)
    tens = capstan_tensions(bent, 100.0)
    check("tension starts at the applied force", tens[1], 100.0)
    check("tension never exceeds the applied force", maximum(tens) <= 100.0, true)
    check("tension never goes negative", minimum(tens) >= 0.0, true)
    check("tension is non-increasing along the rod",
          all(tens[i + 1] <= tens[i] + 1e-12 for i in 1:(length(tens) - 1)), true)
    check("a curved rod loses tension in the core", minimum(tens) < 100.0, true)

    println("\n[4] the breaking-force solve is exact")
    # safety_factor(rod, breaking_force(rod)) must be exactly 1: the load at
    # which some sample first reaches UTS.
    for (n, R, d) in ((60, 10.0, 1.0), (80, 10.0, 1.0), (60, 12.0, 1.0))
        rod = Rod(trefoil_rod(n, R); diameter = d)
        bf = breaking_force(rod)
        check("trefoil R=$R d=$d: breaking force is positive", bf > 0, true)
        check_close("trefoil R=$R d=$d: safety factor at the breaking force",
                    safety_factor(rod, bf), 1.0; atol = 1e-9)
    end

    println("\n[5] knot efficiency")
    # A straight rod is as strong as the material: efficiency exactly 1.
    st1 = Rod(straight(6); diameter = 1.0)
    check_close("straight rod efficiency", knot_efficiency(st1), 1.0; atol = 1e-12)
    check_close("straight rod breaking force = uts * area",
                breaking_force(st1), st1.uts_mpa * pi * (st1.diameter / 2)^2; rtol = 1e-12)
    # A knotted rod is weaker, and lands in the documented 0.60-0.75 band.
    knotted = Rod(trefoil_rod(60, 10.0); diameter = 1.0, friction = 0.35)
    eff = knot_efficiency(knotted)
    check("knotted rod is weaker than the straight strand", eff < 1.0, true)
    check("knotted rod efficiency is positive", eff > 0.0, true)
    check("knotted rod efficiency is in the documented 0.60-0.75 band",
          0.60 < eff < 0.75, true)
    # A strand too thick for the knot's curvature fails from bending alone.
    thick = Rod(trefoil_rod(60, 10.0); diameter = 4.0)
    check("too-thick strand breaks from bending alone (zero tensile capacity)",
          breaking_force(thick), 0.0)
    check_close("its safety factor at zero load is uts / peak bending stress",
                safety_factor(thick, 0.0),
                thick.uts_mpa / maximum(stress_profile(thick, 0.0).total); atol = 1e-9)

    println("\n[6] writhe and linking (corrected quadrature)")
    planar = circle(64, 5.0)
    check_close("planar closed curve has zero writhe", writhe(planar), 0.0; atol = 1e-12)
    check_close("a short curve has zero writhe", writhe([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]), 0.0; atol = 0.0)
    a, b = hopf_pair(64)
    check_close("Hopf pair links once", linking(a, b), -1.0; atol = 0.01)
    check_close("linking is symmetric in its arguments",
                linking(b, a), linking(a, b); atol = 1e-9)
    check_close("reversing one component negates the linking number",
                linking(a, reverse(b)), -linking(a, b); atol = 0.01)
    check_close("reversing both leaves it unchanged",
                linking(reverse(a), reverse(b)), linking(a, b); atol = 1e-9)
    check_close("writhe is orientation-independent",
                writhe(reverse(planar)), writhe(planar); atol = 1e-12)
    check_close("trefoil centreline writhe is close to 3",
                writhe(trefoil_rod(120, 10.0)), 3.0; atol = 0.5)

    println("\n[7] the defect in src/mechanics.jl, pinned")
    # The unnormalised form converges to zero: a genuine Hopf pair reads as
    # unlinked, and the writhe is not even stable under refinement.
    check("unnormalised quadrature misses the Hopf link",
          abs(linking_unnormalised(a, b)) < 0.1, true)
    check("while the corrected form finds it",
          abs(linking(a, b) + 1.0) < 0.01, true)
    unstable = [writhe_unnormalised(trefoil_rod(n, 10.0)) for n in (80, 120, 240)]
    check("unnormalised trefoil writhe is not stable under refinement",
          maximum(abs.(unstable .- writhe(trefoil_rod(120, 10.0)))) > 0.5, true)
    stable = [writhe(trefoil_rod(n, 10.0)) for n in (80, 120, 240)]
    check("corrected trefoil writhe is stable under refinement",
          maximum(abs.(stable .- writhe(trefoil_rod(240, 10.0)))) < 0.1, true)

    println("\n[8] Calugareanu-White-Fuller")
    c0 = calugareanu_check(trefoil_rod(120, 10.0); twist_turns = 0.0)
    check_close("zero-twist residual is small", c0.residual, 0.0; atol = 0.05)
    check_close("zero-twist linking equals the writhe", c0.linking, c0.writhe; atol = 0.05)
    c1 = calugareanu_check(trefoil_rod(120, 10.0); twist_turns = 1.0)
    check("the radial frame is not a material frame: nonzero twist does not close",
          abs(c1.residual) > 0.5, true)

    println()
    if FAILURES[] == 0
        println("ALL PASS — the mechanics-kernel prototype is sound.")
    else
        println("$(FAILURES[]) FAILURE(S) — the mechanics-kernel prototype is wrong.")
        exit(1)
    end
end

main()

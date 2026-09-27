# SPDX-License-Identifier: MPL-2.0
# protistology.jl — minimal working prototype for the protistology kernel
# (src/protistology.jl).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/protistology.jl`; it exits
# non-zero unless every pinned assertion holds.
#
# What this proves: the three flagellar waveforms have the curvature
# statistics they should, and the kinetoplast catenane really is a chainmail
# chain — adjacent rings link exactly once and next-nearest rings link zero,
# but only while the spacing sits inside the window `r < spacing < 2r`.
#
# That last part depends on the Gauss linking integral, and it is where the
# `linking_of` defect documented on feature/knot-mechanics bites hardest:
# with the unnormalised tangents in src/mechanics.jl the whole linking matrix
# comes out as zero, so the chainmail condition is invisible. The prototype
# uses the corrected quadrature (unit tangents, closed-polygon arc lengths)
# and keeps `linking_unnormalised` to demonstrate the difference.

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

segment_lengths(p) = [_norm(_sub(p[i + 1], p[i])) for i in 1:(length(p) - 1)]
segment_lengths_closed(p) =
    [_norm(_sub(p[mod1(i + 1, length(p))], p[i])) for i in eachindex(p)]
total_length(p) = sum(segment_lengths(p))

"""Discrete curvature `|v x a| / |v|^3`, end samples mirroring their neighbours."""
function curvature_profile(p)
    n = length(p)
    n >= 3 || error("curvature_profile needs at least three samples")
    kappa = zeros(Float64, n)
    for i in 2:(n - 1)
        v = _scale(0.5, _sub(p[i + 1], p[i - 1]))
        a = _sub(_add(p[i + 1], p[i - 1]), _scale(2.0, p[i]))
        vlen = _norm(v)
        kappa[i] = vlen > 1e-12 ? _norm(_cross(v, a)) / vlen^3 : 0.0
    end
    kappa[1] = kappa[2]
    kappa[n] = kappa[n - 1]
    kappa
end

# --- flagella ---------------------------------------------------------------

"""
    flagellum(kind; length, samples, amplitude, wavelength, phase, radius)

Centreline of a model flagellum in micrometres. `:planar` is a sinusoidal
planar wave, `:helical` a helix of the given radius and wavelength, `:tipped`
a planar wave with an amplitude envelope growing from base to tip.
"""
function flagellum(kind::Symbol; length::Float64 = 50.0, samples::Int = 200,
                   amplitude::Float64 = 4.0, wavelength::Float64 = 20.0,
                   phase::Float64 = 0.0, radius::Float64 = 2.0)
    samples >= 3 || error("flagellum: need at least 3 samples")
    wavelength > 0 || error("flagellum: wavelength must be positive")
    pts = Vector{NTuple{3,Float64}}(undef, samples)
    for i in 1:samples
        s = length * (i - 1) / (samples - 1)
        if kind == :planar
            y = amplitude * sin(2.0 * pi * s / wavelength + phase)
            pts[i] = (s, y, 0.0)
        elseif kind == :helical
            ang = 2.0 * pi * s / wavelength + phase
            pts[i] = (s, radius * cos(ang), radius * sin(ang))
        elseif kind == :tipped
            env = s / length
            y = amplitude * env * sin(2.0 * pi * s / wavelength + phase)
            pts[i] = (s, y, 0.0)
        else
            error("flagellum: unknown kind :$kind (use :planar, :helical, :tipped)")
        end
    end
    pts
end

"""Geometric summary of one beat frame."""
function beat_summary(points::Vector{NTuple{3,Float64}})
    kappa = curvature_profile(points)
    segs = segment_lengths(points)
    total_k = 0.0
    bend_load = 0.0
    for i in eachindex(segs)
        kbar = 0.5 * (kappa[i] + kappa[i + 1])
        total_k += kbar * segs[i]
        bend_load += kbar^2 * segs[i]
    end
    (samples = length(points), arclength = total_length(points),
     kappa_max = maximum(kappa), kappa_mean = sum(kappa) / length(kappa),
     total_curvature = total_k, bending_load = bend_load)
end

# --- kinetoplast catenanes --------------------------------------------------

"""
    kinetoplast_chain(count; circle_radius, spacing, samples)

A catenane chain: circle `i` lies in the xy-plane centred at
`((i-1) * spacing, 0, 0)` and circle `i+1` threads it in the xz-plane. The
spacing window `r < spacing < 2r` is the chainmail condition — adjacent rings
link once, next-nearest rings link zero.
"""
function kinetoplast_chain(count::Int; circle_radius::Float64 = 1.0,
                           spacing::Float64 = 1.5, samples::Int = 96)
    count >= 1 || error("kinetoplast_chain: count must be >= 1")
    chain = Vector{Vector{NTuple{3,Float64}}}(undef, count)
    for i in 1:count
        ring = Vector{NTuple{3,Float64}}(undef, samples)
        cx = (i - 1) * spacing
        for j in 1:samples
            ang = 2.0 * pi * (j - 1) / samples
            if isodd(i)
                ring[j] = (cx + circle_radius * cos(ang), circle_radius * sin(ang), 0.0)
            else
                ring[j] = (cx + circle_radius * cos(ang), 0.0, circle_radius * sin(ang))
            end
        end
        chain[i] = ring
    end
    chain
end

# --- linking (corrected quadrature) ----------------------------------------

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

"""The src/mechanics.jl form: raw central differences, open-polygon length."""
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

"""Pairwise linking numbers, symmetric with a zero diagonal."""
function linking_matrix(rings::Vector{Vector{NTuple{3,Float64}}})
    n = length(rings)
    m = zeros(Float64, n, n)
    for i in 1:n, j in (i + 1):n
        lk = linking(rings[i], rings[j])
        m[i, j] = lk
        m[j, i] = lk
    end
    m
end

# --- the pins ---------------------------------------------------------------

function main()
    println("protistology prototype")
    println("  (flagellar waveforms, beat summaries, kinetoplast chainmail)")

    println("\n[1] flagellar waveforms validate")
    ok = false
    try
        flagellum(:bogus)
    catch
        ok = true
    end
    check("an unknown kind raises", ok, true)
    ok = false
    try
        flagellum(:planar; samples = 2)
    catch
        ok = true
    end
    check("fewer than three samples raises", ok, true)
    ok = false
    try
        flagellum(:planar; wavelength = 0.0)
    catch
        ok = true
    end
    check("a non-positive wavelength raises", ok, true)

    println("\n[2] planar wave: kappa_max = A k^2")
    # y = A sin(ks): the curvature peaks where the slope vanishes, so
    # kappa_max = A k^2 exactly (up to the finite-difference error).
    k = 2pi / 20.0
    planar = beat_summary(flagellum(:planar))
    check_close("planar kappa_max ~ A k^2", planar.kappa_max, 4.0 * k^2; rtol = 1e-3)
    check("planar wave is longer than its axis (the sine stretches it)",
          planar.arclength > 50.0, true)
    check_close("planar arclength", planar.arclength, 66.0259; rtol = 1e-3)
    check_close("planar kappa_mean", planar.kappa_mean, 0.155806; rtol = 1e-3)
    check_close("planar total curvature", planar.total_curvature, 8.9949; rtol = 1e-3)
    check_close("planar bending load", planar.bending_load, 2.2162; rtol = 1e-3)
    check("planar samples", planar.samples, 200)

    println("\n[3] helical wave: curvature is constant, and equals r k^2/(1+(rk)^2)")
    helical = beat_summary(flagellum(:helical))
    kap = curvature_profile(flagellum(:helical))
    interior = view(kap, 5:196)
    check_close("helical curvature is constant along the interior",
                maximum(interior) - minimum(interior), 0.0; atol = 1e-9)
    check_close("helical kappa = r k^2 / (1 + (r k)^2)",
                helical.kappa_max, 2.0 * k^2 / (1.0 + (2.0 * k)^2); rtol = 1e-4)
    check_close("helical kappa_mean equals kappa_max (constant curvature)",
                helical.kappa_mean, helical.kappa_max; rtol = 1e-12)
    check_close("helical arclength = length * sqrt(1 + (r k)^2)",
                helical.arclength, 50.0 * sqrt(1.0 + (2.0 * k)^2); rtol = 1e-3)
    check_close("helical total curvature", helical.total_curvature, 8.3569; rtol = 1e-3)
    check_close("helical bending load", helical.bending_load, 1.1828; rtol = 1e-3)
    check("the helix bends less than the planar wave",
          helical.bending_load < planar.bending_load, true)

    println("\n[4] tipped wave: envelope grows, so it peaks below the planar wave")
    tipped = beat_summary(flagellum(:tipped))
    check_close("tipped kappa_max", tipped.kappa_max, 0.359569; rtol = 1e-3)
    check("tipped peaks below the planar wave", tipped.kappa_max < planar.kappa_max, true)
    check("tipped arclength is shorter than the planar one",
          tipped.arclength < planar.arclength, true)
    check_close("tipped arclength", tipped.arclength, 55.7839; rtol = 1e-3)
    check_close("tipped bending load", tipped.bending_load, 0.9433; rtol = 1e-3)
    # the phase shift must move the waveform without changing its statistics
    shifted = beat_summary(flagellum(:planar; phase = 1.0))
    check_close("phase does not change kappa_max",
                shifted.kappa_max, planar.kappa_max; rtol = 1e-9)
    check_close("phase does not change the arclength",
                shifted.arclength, planar.arclength; rtol = 1e-9)

    println("\n[5] kinetoplast chain: the chainmail condition")
    rings = kinetoplast_chain(4)
    M = linking_matrix(rings)
    check("a single ring links nothing", linking_matrix(kinetoplast_chain(1)), zeros(1, 1))
    check_close("adjacent rings link once (1,2)", M[1, 2], -1.0; atol = 0.01)
    check_close("adjacent rings link once (2,3)", M[2, 3], 1.0; atol = 0.01)
    check_close("adjacent rings link once (3,4)", M[3, 4], -1.0; atol = 0.01)
    check_close("next-nearest rings link zero (1,3)", M[1, 3], 0.0; atol = 0.01)
    check_close("next-nearest rings link zero (2,4)", M[2, 4], 0.0; atol = 0.01)
    check_close("far rings link zero (1,4)", M[1, 4], 0.0; atol = 0.01)
    check("the matrix is symmetric", M, transpose(M))
    check_close("the diagonal is zero", sum(abs(M[i, i]) for i in axes(M, 1)), 0.0; atol = 1e-12)
    check("the linking signs alternate along the chain",
          sign(M[1, 2]) != sign(M[2, 3]) && sign(M[2, 3]) != sign(M[3, 4]), true)
    # five rings: still a path, still alternating
    M5 = linking_matrix(kinetoplast_chain(5))
    check_close("five rings: (1,2) links", M5[1, 2], -1.0; atol = 0.01)
    check_close("five rings: (4,5) links", M5[4, 5], 1.0; atol = 0.01)
    check_close("five rings: (1,3) does not", M5[1, 3], 0.0; atol = 0.01)
    check_close("five rings: (1,5) does not", M5[1, 5], 0.0; atol = 0.01)

    println("\n[6] the chainmail window r < spacing < 2r")
    # Inside the window adjacent rings thread; outside it they miss.
    for spacing in (1.2, 1.5, 1.8)
        r2 = kinetoplast_chain(2; spacing = spacing)
        check_close("spacing $spacing (inside the window): rings link",
                    abs(linking(r2[1], r2[2])), 1.0; atol = 0.01)
    end
    for spacing in (2.5, 3.0)
        r2 = kinetoplast_chain(2; spacing = spacing)
        check_close("spacing $spacing (outside the window): rings miss",
                    abs(linking(r2[1], r2[2])), 0.0; atol = 0.01)
    end
    ok = false
    try
        kinetoplast_chain(0)
    catch
        ok = true
    end
    check("a chain of zero rings raises", ok, true)

    println("\n[7] what the linking defect costs")
    # With the unnormalised tangents in src/mechanics.jl the entire matrix is
    # zero, so the chainmail condition is invisible.
    rings4 = kinetoplast_chain(4)
    defective = maximum(abs, [linking_unnormalised(rings4[i], rings4[j])
                              for i in 1:4, j in 1:4 if i != j])
    corrected = maximum(abs, [linking(rings4[i], rings4[j])
                              for i in 1:4, j in 1:4 if i != j])
    check("the corrected quadrature finds the links", corrected > 0.9, true)
    check("the unnormalised quadrature loses them entirely", defective < 0.1, true)

    println()
    if FAILURES[] == 0
        println("ALL PASS — the protistology-kernel prototype is sound.")
    else
        println("$(FAILURES[]) FAILURE(S) — the protistology-kernel prototype is wrong.")
        exit(1)
    end
end

main()

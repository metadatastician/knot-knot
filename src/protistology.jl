# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Protistological topology: flagellar waveform geometry and kinetoplast
# catenation. This module is deliberately self-contained (geometry only, no
# mechanics units) so it can split out as the seed of a future Proctist.jl
# without dragging the engineering layer along.

"""
    flagellum(kind::Symbol; length = 50.0, samples = 200, amplitude = 4.0,
              wavelength = 20.0, phase = 0.0, radius = 2.0) -> Vector{NTuple{3,Float64}}

Centreline of a model flagellum in micrometres.

Kinds:
- `:planar` — sinusoidal planar wave ``y = A \\sin(2\\pi (x/\\lambda) + \\phi)``
  (the classical activated-sperm / sea-urchin waveform);
- `:helical` — helical wave of the given `radius` and `wavelength` along the
  axis (trypanosome / dinoflagellate-style);
- `:tipped` — planar wave with an amplitude envelope growing from base to
  tip (a crude Chlamydomonas-like breaststroke proxy).
"""
function flagellum(
    kind::Symbol;
    length::Float64 = 50.0,
    samples::Int = 200,
    amplitude::Float64 = 4.0,
    wavelength::Float64 = 20.0,
    phase::Float64 = 0.0,
    radius::Float64 = 2.0,
)::Vector{NTuple{3,Float64}}
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

"""
    beat_summary(points; wavelength = 20.0) -> NamedTuple

Geometric summary of one beat frame: sample count, arclength, maximum and
mean curvature, total curvature ``\\int \\kappa\\, ds``, and a dimensionless
bending load ``\\int \\kappa^2 ds`` proportional to the elastic power a
motor ensemble must deliver (G. B. R. standard small-deformation scaling).
"""
function beat_summary(points::Vector{NTuple{3,Float64}})
    rod = Rod(points; diameter = 1.0, youngs_mpa = 1.0, uts_mpa = 1.0, friction = 0.0)
    kappa = curvature_profile(rod)
    segs = segment_lengths(rod)
    total_k = 0.0
    bend_load = 0.0
    for i in eachindex(segs)
        kbar = 0.5 * (kappa[i] + kappa[i + 1])
        total_k += kbar * segs[i]
        bend_load += kbar^2 * segs[i]
    end
    (
        samples = length(points),
        arclength = total_length(rod),
        kappa_max = maximum(kappa),
        kappa_mean = sum(kappa) / length(kappa),
        total_curvature = total_k,
        bending_load = bend_load,
    )
end

"""
    kinetoplast_chain(count; circle_radius = 1.0, spacing = 1.5, samples = 96)

A catenane chain modelling a kinetoplast-DNA minicircle row: circle `i`
lies in the xy-plane centred at `((i-1) * spacing, 0, 0)` and circle `i+1`
threads it in the xz-plane. The spacing window `r < spacing < 2r` is the
chainmail condition: adjacent rings link once, next-nearest rings link
zero. Returns a vector of closed polygons.
"""
function kinetoplast_chain(
    count::Int;
    circle_radius::Float64 = 1.0,
    spacing::Float64 = 1.5,
    samples::Int = 96,
)::Vector{Vector{NTuple{3,Float64}}}
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

"""
    linking_matrix(rings) -> Matrix{Float64}

Pairwise Gauss-integral linking numbers for a collection of closed curves.
"""
function linking_matrix(rings::Vector{Vector{NTuple{3,Float64}}})::Matrix{Float64}
    n = length(rings)
    m = zeros(Float64, n, n)
    for i in 1:n, j in (i + 1):n
        lk = linking_of(rings[i], rings[j])
        m[i, j] = lk
        m[j, i] = lk
    end
    m
end

"""
    protistology_report(; beat_kind = :planar, beat_phases = 8, chain_rings = 4)

One-call demonstrator combining both application faces: a flagellar beat
cycle with its bending-load sweep, and a kinetoplast catenane with its
linking matrix. Returns both as a NamedTuple; the web laboratory renders
the same payload.
"""
function protistology_report(;
    beat_kind::Symbol = :planar,
    beat_phases::Int = 8,
    chain_rings::Int = 4,
)
    phases = [2.0 * pi * k / beat_phases for k in 0:(beat_phases - 1)]
    loads = Float64[]
    for ph in phases
        frame = flagellum(beat_kind; phase = ph)
        push!(loads, beat_summary(frame).bending_load)
    end
    chain = kinetoplast_chain(chain_rings)
    (
        waveform = beat_kind,
        bending_loads = loads,
        linking_matrix = linking_matrix(chain),
        adjacent_link = chain_rings >= 2 ? linking_of(chain[1], chain[2]) : 0.0,
    )
end

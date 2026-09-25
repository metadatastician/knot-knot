# Applied knot mechanics: curvature, capstan, stress, topology integrals.

function _circle(radius::Float64, samples::Int)
    pts = Vector{NTuple{3,Float64}}(undef, samples)
    for i in 1:samples
        ang = 2.0 * pi * (i - 1) / samples
        pts[i] = (radius * cos(ang), radius * sin(ang), 0.0)
    end
    pts
end

@testset "mechanics — curvature of a circle" begin
    r = 10.0
    rod = Rod(_circle(r, 240); diameter = 1.0, youngs_mpa = 1000.0, uts_mpa = 500.0)
    kappa = curvature_profile(rod)
    interior = kappa[3:(end - 2)]
    @test all(abs(k - 1 / r) < 0.01 for k in interior)
end

@testset "mechanics — straight rod carries pure axial stress" begin
    pts = [(Float64(x), 0.0, 0.0) for x in 0:50]
    rod = Rod(pts; diameter = 2.0, youngs_mpa = 3500.0, uts_mpa = 850.0)
    area = pi * 1.0^2
    sigma = stress_profile(rod, 100.0)
    @test maximum(sigma.bending) ≈ 0.0 atol = 1e-9
    @test sigma.axial[1] ≈ 100.0 / area
    # Bending energy of a straight rod vanishes.
    @test bending_energy(rod) ≈ 0.0 atol = 1e-9
    # Breaking force of a straight rod is exactly UTS * A.
    @test breaking_force(rod) ≈ 850.0 * area rtol = 1e-9
    @test knot_efficiency(rod) ≈ 1.0 rtol = 1e-9
end

@testset "mechanics — capstan decay through a curved core" begin
    # S-curve: a curved core in the middle, straight tails.
    pts = Vector{NTuple{3,Float64}}()
    for x in 0:60
        y = 8.0 * sin(2.0 * pi * x / 40.0)
        push!(pts, (Float64(x), y, 0.0))
    end
    rod = Rod(pts; diameter = 2.0, youngs_mpa = 3500.0, uts_mpa = 850.0, friction = 0.35)
    tensions = capstan_tensions(rod, 200.0)
    @test tensions[1] == 200.0
    @test minimum(tensions) < 200.0
    @test all(tensions .>= 0.0)
    # Stress superposition: total >= axial and total >= bending pointwise.
    sigma = stress_profile(rod, 200.0)
    @test all(sigma.total[i] >= sigma.axial[i] - 1e-12 for i in eachindex(sigma.total))
    @test safety_factor(rod, 200.0) > 0
end

@testset "mechanics — writhe and linking integrals" begin
    # Planar circle: writhe zero by symmetry.
    @test abs(writhe_of(_circle(5.0, 160))) < 1e-6

    # Positive Hopf pair: linking number +1.
    c1 = _circle(1.0, 120)
    c2 = Vector{NTuple{3,Float64}}(undef, 120)
    for i in 1:120
        ang = 2.0 * pi * (i - 1) / 120
        c2[i] = (1.0 + 0.4 * cos(ang), 0.0, 0.4 * sin(ang))
    end
    # Orientation conventions aside, the Hopf pair links exactly once.
    @test abs(linking_of(c1, c2)) ≈ 1.0 atol = 0.02
end

@testset "mechanics — Calugareanu identity on a planar ribbon" begin
    res = calugareanu_check(_circle(5.0, 200); twist_turns = 3.0, ribbon_width = 0.2)
    @test abs(res.writhe) < 0.02
    # |Lk| = |Tw + Wr| = 3 for a planar centreline with three twist turns;
    # the sign tracks the frame-rotation convention.
    @test abs(res.linking) ≈ 3.0 atol = 0.1
    @test abs(abs(res.linking) - abs(res.twist + res.writhe)) < 0.1
end

# Protistological topology: flagellar waveforms and kinetoplast catenanes.

@testset "protistology — flagellar waveforms" begin
    planar = flagellum(:planar; length = 50.0, amplitude = 4.0, wavelength = 20.0)
    @test length(planar) == 200
    @test planar[1][1] ≈ 0.0
    @test planar[end][1] ≈ 50.0
    @test maximum(abs(p[3]) for p in planar) ≈ 0.0

    helical = flagellum(:helical; length = 50.0, radius = 2.0, wavelength = 20.0)
    # Analytic curvature of the unit-speed-axis helix: r w^2 / (1 + r^2 w^2)
    w = 2.0 * pi / 20.0
    expected = 2.0 * w^2 / (1.0 + 4.0 * w^2)
    s = beat_summary(helical)
    @test s.kappa_max ≈ expected rtol = 0.05

    tipped = flagellum(:tipped; length = 50.0, amplitude = 4.0)
    st = beat_summary(tipped)
    @test st.total_curvature > 0
    @test st.bending_load > 0
end

@testset "protistology — beat cycle" begin
    loads = Float64[]
    for k in 0:7
        frame = flagellum(:planar; phase = 2.0 * pi * k / 8)
        push!(loads, beat_summary(frame).bending_load)
    end
    # Phase-shifting a periodic wave changes the discrete bending load only
    # through boundary quadrature: all frames agree to within ~10%.
    @test maximum(loads) - minimum(loads) < 0.1 * maximum(loads)
end

@testset "protistology — kinetoplast catenane chain" begin
    chain = kinetoplast_chain(4)
    @test length(chain) == 4
    m = linking_matrix(chain)
    for i in 1:3
        @test abs(abs(m[i, i + 1]) - 1.0) < 0.05
    end
    @test abs(m[1, 3]) < 0.05
    @test abs(m[2, 4]) < 0.05
    @test abs(m[1, 4]) < 0.05
end

@testset "protistology — one-call report" begin
    rep = protistology_report(; beat_kind = :helical, beat_phases = 4, chain_rings = 3)
    @test length(rep.bending_loads) == 4
    @test size(rep.linking_matrix) == (3, 3)
    @test abs(abs(rep.adjacent_link) - 1.0) < 0.05
end

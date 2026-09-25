# Diagram kernel: PD conventions, writhe, mirror, Gauss traversals.

@testset "diagrams — kernel conventions" begin
    tre = standard_knot("3_1")
    fe = standard_knot("4_1")
    un = standard_knot("0_1")

    @test crossing_count(tre.pd) == 3
    @test writhe(tre.pd) == 3          # right-hand trefoil: 3 positive crossings
    @test writhe(fe.pd) == 0           # figure-eight: alternating, writhe 0
    @test writhe(mirror(tre.pd)) == -3
    @test crossing_count(un.pd) == 0

    @test length(arcs_of(tre.pd)) == 6
    comps = traverse_components(tre.pd)
    @test length(comps) == 1
    @test length(comps[1]) == 6

    # Skein-format Gauss code: every crossing appears once + and once -.
    code = skein_gauss_code(tre.pd)
    @test length(code) == 6
    for k in 1:3
        @test count(x -> x == k, code) == 1
        @test count(x -> x == -k, code) == 1
    end

    # Gauss word visits match crossing signs.
    gw = gauss_word(tre.pd)
    @test length(gw.visits) == 6
    @test all(v.sign == 1 for v in gw.visits)
    @test count(v -> v.over, gw.visits) == 3
end

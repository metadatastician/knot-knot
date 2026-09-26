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

@testset "diagrams — mirror swaps the passages, not just the signs" begin
    # A mirror reflects the projection plane, so at every crossing the two
    # passages that were over become under and vice versa: X[a,b,c,d,s] ->
    # X[a,d,c,b,-s]. Sign negation alone (the pre-#mirror-fix behaviour)
    # leaves a diagram of the wrong knot.
    tre = standard_knot("3_1").pd

    m = mirror(tre)
    @test [c.arcs for c in m.crossings] == [(1, 5, 2, 4), (3, 1, 4, 6), (5, 3, 6, 2)]
    @test all(c.sign == -1 for c in m.crossings)

    # Arc labels survive, so the mirrored diagram is still a valid PD.
    @test Set(arcs_of(m)) == Set(arcs_of(tre))
    @test crossing_count(m) == crossing_count(tre)
    @test writhe(m) == -writhe(tre)

    # Involutive: mirror(mirror(d)) is d again, crossings and signs alike.
    mm = mirror(m)
    @test [c.arcs for c in mm.crossings] == [c.arcs for c in tre.crossings]
    @test [c.sign for c in mm.crossings] == [c.sign for c in tre.crossings]

    # The unknot has no crossings, so it mirrors to itself.
    @test mirror(standard_knot("0_1").pd).crossings == Crossing[]
end

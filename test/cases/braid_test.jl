# Braid words, closures, and crossing-level sanity.

@testset "braid — parsing and closure" begin
    b = parse_braid_word("s1.s1.s1")
    @test b.generators == [(1, 1), (1, 1), (1, 1)]
    d = braid_closure(b)
    @test crossing_count(d) == 3
    @test writhe(d) == 3

    # The inverse braid is the mirror: signs flip, determinant is unchanged.
    dm = braid_closure(parse_braid_word("S1.S1.S1"))
    @test writhe(dm) == -3
    @test determinant(dm) == 3

    # Figure-eight braid: s1 S2 s1 S2.
    d4 = braid_closure(parse_braid_word("s1.S2.s1.S2"))
    @test crossing_count(d4) == 4
    @test writhe(d4) == 0
    @test genus(d4) == 1

    # Empty braid closes to the unknot.
    @test crossing_count(braid_closure(parse_braid_word(""))) == 0
end

@testset "braid — torus knot T(2,5)" begin
    d = braid_closure(parse_braid_word("s1.s1.s1.s1.s1"))
    @test crossing_count(d) == 5
    @test determinant(d) == 5
    @test genus(d) == 2
end

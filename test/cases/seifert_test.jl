# Seifert circles, genus and signature pins.

@testset "seifert — circles and genus" begin
    @test seifert_circles(standard_knot("3_1").pd) == 2
    @test seifert_circles(standard_knot("4_1").pd) == 3
    @test genus(standard_knot("3_1").pd) == 1
    @test genus(standard_knot("4_1").pd) == 1
    @test genus(standard_knot("5_1").pd) == 2
    @test genus(standard_knot("0_1").pd) == 0
    # KnotAtlas 3-genus values for the rest of the table. The Atlas table
    # diagrams are alternating, so Seifert's algorithm is exact on them.
    for (name, g) in (
        ("5_2", 1),
        ("6_1", 1),
        ("6_2", 2),
        ("6_3", 2),
        ("7_1", 3),
        ("7_2", 1),
        ("7_3", 2),
        ("7_4", 1),
        ("7_5", 2),
        ("7_6", 2),
        ("7_7", 2),
    )
        @test genus(standard_knot(name).pd) == g
    end
    # Genus identity g = (c - s + 1)/2 on every table knot.
    for name in knot_names()
        d = standard_knot(name).pd
        crossing_count(d) == 0 && continue
        @test 2 * genus(d) == crossing_count(d) - seifert_circles(d) + 1
    end
end

@testset "seifert — signature pins" begin
    @test signature(standard_knot("3_1").pd) == -2
    @test signature(standard_knot("4_1").pd) == 0
    @test signature(standard_knot("0_1").pd) == 0
    # Mirror antisymmetry.
    @test signature(mirror(standard_knot("3_1").pd)) == 2
end

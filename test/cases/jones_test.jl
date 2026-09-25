# Jones polynomial via the Kauffman bracket: pins, mirror, normalisation.

@testset "jones — pins and identities" begin
    tre = standard_knot("3_1").pd
    fe = standard_knot("4_1").pd

    # Right-hand trefoil pin (KnotTheory.jl convention).
    @test jones_polynomial(tre) == LaurentPoly(-4 => -1, -3 => 1, -1 => 1)

    # Unknot normalisation.
    @test jones_polynomial(standard_knot("0_1").pd) == LaurentPoly(0 => 1)

    # V(1) = 1 for every knot.
    for name in knot_names()
        d = standard_knot(name).pd
        v = jones_polynomial(d)
        @test sum(values(v)) == 1
    end

    # Figure-eight is amphichiral: V(t) = V(t^-1).
    vfe = jones_polynomial(fe)
    @test Set(keys(vfe)) == Set(-e for e in keys(vfe))
    for e in keys(vfe)
        @test vfe[e] == vfe[-e]
    end
end

@testset "jones — mirror relation" begin
    tre = standard_knot("3_1").pd
    vm = jones_polynomial(mirror(tre))
    vp = jones_polynomial(tre)
    @test Set(keys(vm)) == Set(-e for e in keys(vp))
    for e in keys(vp)
        @test vm[-e] == vp[e]
    end
end

@testset "jones — KnotAtlas pins across the table" begin
    # Published Jones polynomials (KnotAtlas) for every nontrivial table
    # knot. These pin both the sign assignments of the table PDs and the
    # bracket normalisation V = (-A^3)^w <D> at A = t^(-1/4).
    atlas = Dict(
        "5_1" => LaurentPoly(-7 => -1, -6 => 1, -5 => -1, -4 => 1, -2 => 1),
        "5_2" => LaurentPoly(-6 => -1, -5 => 1, -4 => -1, -3 => 2, -2 => -1, -1 => 1),
        "6_1" => LaurentPoly(-4 => 1, -3 => -1, -2 => 1, -1 => -2, 0 => 2, 1 => -1, 2 => 1),
        "6_2" => LaurentPoly(
            -5 => 1,
            -4 => -2,
            -3 => 2,
            -2 => -2,
            -1 => 2,
            0 => -1,
            1 => 1,
        ),
        "6_3" => LaurentPoly(-3 => -1, -2 => 2, -1 => -2, 0 => 3, 1 => -2, 2 => 2, 3 => -1),
        "7_1" => LaurentPoly(
            -10 => -1,
            -9 => 1,
            -8 => -1,
            -7 => 1,
            -6 => -1,
            -5 => 1,
            -3 => 1,
        ),
        "7_2" => LaurentPoly(
            -8 => -1,
            -7 => 1,
            -6 => -1,
            -5 => 2,
            -4 => -2,
            -3 => 2,
            -2 => -1,
            -1 => 1,
        ),
        "7_3" => LaurentPoly(
            2 => 1,
            3 => -1,
            4 => 2,
            5 => -2,
            6 => 3,
            7 => -2,
            8 => 1,
            9 => -1,
        ),
        "7_4" => LaurentPoly(
            1 => 1,
            2 => -2,
            3 => 3,
            4 => -2,
            5 => 3,
            6 => -2,
            7 => 1,
            8 => -1,
        ),
        "7_5" => LaurentPoly(
            -9 => -1,
            -8 => 2,
            -7 => -3,
            -6 => 3,
            -5 => -3,
            -4 => 3,
            -3 => -1,
            -2 => 1,
        ),
        "7_6" => LaurentPoly(
            -9 => -1,
            -8 => 2,
            -7 => -3,
            -6 => 4,
            -5 => -3,
            -4 => 3,
            -3 => -2,
            -2 => 1,
        ),
        "7_7" => LaurentPoly(
            -3 => -1,
            -2 => 3,
            -1 => -3,
            0 => 4,
            1 => -4,
            2 => 3,
            3 => -2,
            4 => 1,
        ),
    )
    for (name, want) in atlas
        @test jones_polynomial(standard_knot(name).pd) == want
    end
end

@testset "jones — braid closures" begin
    trefoil_braid = braid_closure(parse_braid_word("s1.s1.s1"))
    @test jones_polynomial(trefoil_braid) == LaurentPoly(-4 => -1, -3 => 1, -1 => 1)
    fig8_braid = braid_closure(parse_braid_word("s1.S2.s1.S2"))
    @test determinant(fig8_braid) == 5
    # A closed braid's Jones polynomial always satisfies V(1) = 1.
    @test sum(values(jones_polynomial(fig8_braid))) == 1
end

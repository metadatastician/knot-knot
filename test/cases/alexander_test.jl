# Alexander-Conway pins from the literature (Conway normalisation).

@testset "alexander — pinned polynomials" begin
    tre = standard_knot("3_1").pd
    fe = standard_knot("4_1").pd
    st = standard_knot("6_1").pd

    @test alexander_polynomial(tre) == LaurentPoly(-1 => 1, 0 => -1, 1 => 1)
    @test alexander_polynomial(fe) == LaurentPoly(-1 => -1, 0 => 3, 1 => -1)
    @test alexander_polynomial(standard_knot("0_1").pd) == LaurentPoly(0 => 1)

    # Symmetry Delta(t) = Delta(t^-1) for knots.
    for d in (tre, fe, st)
        ax = alexander_polynomial(d)
        @test Set(keys(ax)) == Set(-e for e in keys(ax))
    end

    # Conway normalisation: Delta(1) = +1.
    @test sum(values(alexander_polynomial(tre))) == 1
    @test sum(values(alexander_polynomial(fe))) == 1
end

@testset "alexander — KnotAtlas pins across the table" begin
    # Published Alexander polynomials (KnotAtlas), Conway normalisation.
    # These pin the orientation-sensitive Fox row: the exponent at each
    # crossing is fixed by the over-strand's entering slot, not by the
    # crossing sign (6_3 mixes both configurations).
    atlas = Dict(
        "5_1" => LaurentPoly(-2 => 1, -1 => -1, 0 => 1, 1 => -1, 2 => 1),
        "5_2" => LaurentPoly(-1 => 2, 0 => -3, 1 => 2),
        "6_1" => LaurentPoly(-1 => -2, 0 => 5, 1 => -2),
        "6_2" => LaurentPoly(-2 => -1, -1 => 3, 0 => -3, 1 => 3, 2 => -1),
        "6_3" => LaurentPoly(-2 => 1, -1 => -3, 0 => 5, 1 => -3, 2 => 1),
        "7_1" => LaurentPoly(-3 => 1, -2 => -1, -1 => 1, 0 => -1, 1 => 1, 2 => -1, 3 => 1),
        "7_2" => LaurentPoly(-1 => 3, 0 => -5, 1 => 3),
        "7_3" => LaurentPoly(-2 => 2, -1 => -3, 0 => 3, 1 => -3, 2 => 2),
        "7_4" => LaurentPoly(-1 => 4, 0 => -7, 1 => 4),
        "7_5" => LaurentPoly(-2 => 2, -1 => -4, 0 => 5, 1 => -4, 2 => 2),
        "7_6" => LaurentPoly(-2 => -1, -1 => 5, 0 => -7, 1 => 5, 2 => -1),
        "7_7" => LaurentPoly(-2 => 1, -1 => -5, 0 => 9, 1 => -5, 2 => 1),
    )
    for (name, want) in atlas
        @test alexander_polynomial(standard_knot(name).pd) == want
    end
end

@testset "alexander — determinants across the table" begin
    # KnotAtlas "Determinant and Signature" values for every table knot.
    known = Dict(
        "0_1" => 1,
        "3_1" => 3,
        "4_1" => 5,
        "5_1" => 5,
        "5_2" => 7,
        "6_1" => 9,
        "6_2" => 11,
        "6_3" => 13,
        "7_1" => 7,
        "7_2" => 11,
        "7_3" => 13,
        "7_4" => 15,
        "7_5" => 17,
        "7_6" => 19,
        "7_7" => 21,
    )
    for (name, det) in known
        @test determinant(standard_knot(name).pd) == det
    end
    # Jones at -1 agrees with |Delta(-1)| on every table knot.
    for name in knot_names()
        d = standard_knot(name).pd
        crossing_count(d) == 0 && continue
        v = jones_polynomial(d)
        @test abs(sum(c * (-1)^e for (e, c) in v)) == determinant(d)
    end
end

@testset "conway — pinned polynomials" begin
    @test conway_polynomial(standard_knot("3_1").pd) == LaurentPoly(0 => 1, 2 => 1)
    @test conway_polynomial(standard_knot("4_1").pd) == LaurentPoly(0 => 1, 2 => -1)
    @test conway_polynomial(standard_knot("0_1").pd) == LaurentPoly(0 => 1)
end

# Skein relations: switching, oriented smoothing, and the Conway identity.

@testset "skein — Conway relation holds on the trefoil" begin
    d = standard_knot("3_1").pd
    for i in 1:3
        @test iszero(verify_conway_skein(d, i))
    end
end

@testset "skein — Conway relation holds on the figure-eight" begin
    d = standard_knot("4_1").pd
    for i in 1:4
        @test iszero(verify_conway_skein(d, i))
    end
end

@testset "skein — trefoil triple values" begin
    d = standard_knot("3_1").pd
    lp, lm, l0 = skein_triple(d, 1)
    # L+ is the trefoil, L- the unknot, L0 the Hopf link.
    @test conway_polynomial(lp) == LaurentPoly(0 => 1, 2 => 1)
    @test conway_polynomial(lm) == LaurentPoly(0 => 1)
    n0 = conway_polynomial(l0)
    @test n0 == LaurentPoly(1 => 1) || n0 == LaurentPoly(1 => -1)
end

@testset "skein — smoothing reduces crossings" begin
    d = standard_knot("3_1").pd
    l0 = smooth_crossing(d, 1)
    @test crossing_count(l0) == 2
    # The smoothing preserves every other crossing's sign.
    @test all(c.sign == 1 for c in l0.crossings)
end

@testset "skein — Hopf link from a braid" begin
    hopf = braid_closure(parse_braid_word("s1.s1"))
    @test crossing_count(hopf) == 2
    n = conway_polynomial(hopf)
    @test n == LaurentPoly(1 => 1) || n == LaurentPoly(1 => -1)
end

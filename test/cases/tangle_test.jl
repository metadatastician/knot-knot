# Conway's rational tangle calculus: fractions, continued fractions,
# 2-bridge classification.

@testset "tangle — continued fractions" begin
    @test continued_fraction([3]) == 3 // 1
    # 1 + 1/(1 + 1/(1 + 1/1)) = 5/3 — the figure-eight fraction.
    @test continued_fraction([1, 1, 1, 1]) == 5 // 3
    @test continued_fraction([2, 2]) == 5 // 2

    # Round trip: notation -> fraction -> notation evaluates identically.
    for terms in ([3], [1, 2], [-2, 3], [1, 1, 1, 1], [2, -1, 4])
        f = continued_fraction(terms)
        @test continued_fraction(conway_notation(f)) == f
    end
end

@testset "tangle — arithmetic" begin
    t3 = RationalTangle([3])
    t1 = RationalTangle([1])
    @test tangle_fraction(t3 + t1) == 4 // 1
    @test tangle_fraction(t3 * t1) == 3 // 1
    @test tangle_fraction(mirror(t3)) == -3 // 1
    @test tangle_fraction(rotate(t1)) == -1 // 1
end

@testset "tangle — two-bridge classification" begin
    @test two_bridge_det(3 // 1) == 3          # trefoil
    @test two_bridge_det(5 // 3) == 5          # figure-eight
    @test two_bridge_det(5 // 2) == 5          # three-twist knot
    @test two_bridge_det(9 // 2) == 9          # stevedore
    @test is_knot_fraction(3 // 1)
    @test is_knot_fraction(5 // 2)
    @test !is_knot_fraction(2 // 1)            # even numerator: 2-component link
    @test !is_knot_fraction(4 // 3)
end

# Laurent polynomial arithmetic — the substrate every invariant builds on.

@testset "laurent — arithmetic" begin
    a = LaurentPoly(0 => 1, 1 => 2)
    b = LaurentPoly(-1 => 1, 1 => -2)
    @test haskey(a + b, 0)
    @test (a + b)[0] == 2
    @test !haskey(a + b, 1)
    @test (a * b)[0] == -4 + 1
    @test lp_pow(LaurentPoly(1 => 1), 3) == LaurentPoly(3 => 1)
    @test lp_shift(a, -2) == LaurentPoly(-2 => 1, -1 => 2)
    @test lp_eval(LaurentPoly(-1 => 1, 1 => 1), 2) == 2 + 1 // 2
    @test lp_show_string(LaurentPoly()) == "0"
    s = lp_show_string(LaurentPoly(-4 => -1, -3 => 1, -1 => 1))
    @test s == "-t^-4 + t^-3 + t^-1"
end

# Quandle axioms and colouring-count pins.

@testset "quandle — axioms" begin
    @test is_quandle(dihedral_quandle(3))
    @test is_quandle(dihedral_quandle(5))
    @test is_quandle(alexander_quandle(4, 3))
    # A deliberately broken table: idempotence fails.
    bad = Quandle([2 1; 1 2])
    @test !is_quandle(bad)
end

@testset "quandle — colouring counts" begin
    r3 = dihedral_quandle(3)
    r5 = dihedral_quandle(5)
    @test coloring_count(standard_knot("0_1").pd, r3) == 3
    @test coloring_count(standard_knot("3_1").pd, r3) == 9
    @test coloring_count(standard_knot("4_1").pd, r3) == 3
    # det(3_1) = 3 is coprime to 5, so only trivial R5 colourings exist.
    @test coloring_count(standard_knot("3_1").pd, r5) == 5
    @test coloring_count(standard_knot("0_1").pd, r5) == 5
end

@testset "quandle — presentation shape" begin
    rels = fundamental_relations(standard_knot("3_1").pd)
    @test length(rels) == 3
    gens = Set{Int}()
    for (x, y, z, s) in rels
        push!(gens, x, y, z)
        @test s == 1
    end
    @test length(gens) == 3
end

# Conway-suite report, JSON bridge, and KnotTheory interop round-trips.

@testset "report — trefoil suite" begin
    r = conway_suite("3_1")
    @test r.crossings == 3
    @test r.writhe == 3
    @test r.determinant == 3
    @test r.signature == -2
    @test r.genus == 1
    @test r.seifert_circles == 2
    @test Dict(r.colorings)["R3"] == 9
    buf = IOBuffer()
    show(buf, r)
    @test occursin("Jones", String(take!(buf)))
end

@testset "report — JSON bridge" begin
    entry = standard_knot("4_1")
    path = tempname() * ".json"
    write_visualization(path, entry)
    text = read(path, String)
    @test occursin("knotknot/visual-v0", text)
    @test occursin("figure-eight", text)
    @test occursin("skein", text)
    @test occursin("quandle", text)
    rm(path)

    rec = export_knot_record(conway_suite("3_1"); gauss = [1, -2, 3, -1, 2, -3])
    @test rec["crossing_number"] == 3
    @test occursin(":", rec["jones_polynomial"])
end

@testset "report — KnotTheory dict round-trip" begin
    d = standard_knot("3_1").pd
    x = to_knottheory_dict(d)
    @test x["convention"] == "knotatlas-pd: under a->c, over d->b, arcs counter-clockwise"
    d2 = from_knottheory_dict(x)
    @test writhe(d2) == writhe(d)
    @test jones_polynomial(d2) == jones_polynomial(d)
    # Skein-format Gauss code is well formed: each label once over, once under.
    code = skein_gauss_code(d)
    for k in 1:3
        @test count(v -> v == k, code) == 1
        @test count(v -> v == -k, code) == 1
    end
end

# SPDX-License-Identifier: MPL-2.0
# export_catalog.jl — write the web-laboratory catalog JSON.
#
# Emits www/public/knot-lab/knot_data.json: one KnotKnot evaluation record
# per standard knot, which the web laboratory layers over its mechanics
# view (quandle, skein and Conway annotations included).

using KnotKnot

outdir = joinpath(@__DIR__, "..", "www", "public", "knot-lab")
mkpath(outdir)

records = Dict{String,Any}[]
for name in knot_names()
    entry = standard_knot(name)
    report = conway_suite(entry)
    push!(
        records,
        Dict{String,Any}(
            "name" => name,
            "record" => export_knot_record(report; gauss = skein_gauss_code(entry.pd)),
            "pd" => [collect(t) for t in pd(entry.pd)],
            "description" => entry.description,
        ),
    )
end

payload = Dict{String,Any}(
    "schema" => "knotknot/catalog-v0",
    "generated_by" => "KnotKnot.jl examples/export_catalog.jl",
    "knots" => records,
    "annotations" => Dict{String,Any}(
        "skein" =>
            "Crossings carry the L+/L-/L0 skein triple; Conway's relation is verified per crossing.",
        "quandle" =>
            "Arcs are fundamental-quandle generators; colouring counts are homomorphism counts.",
        "conway" =>
            "Rational knots are classified by their Conway fraction p/q; determinant = |p|.",
    ),
)

path = joinpath(outdir, "knot_data.json")
open(path, "w") do io
    print(io, to_json_string(payload))
end
println("wrote $(path) with $(length(records)) knots")

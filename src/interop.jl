# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Bridge surfaces: JSON export (stdlib-only serialiser), the Skein.jl
# KnotRecord-compatible interchange shape, and KnotTheory.jl PD-dict
# round-tripping. The serialiser is deliberately small and dependency-free;
# registry JSON packages remain optional downstream choices.

function _json_escape(s::AbstractString)::String
    out = IOBuffer()
    for ch in String(s)
        if ch == '"'
            print(out, "\\\"")
        elseif ch == '\\'
            print(out, "\\\\")
        elseif ch == '\n'
            print(out, "\\n")
        elseif ch == '\t'
            print(out, "\\t")
        elseif ch == '\r'
            print(out, "\\r")
        elseif Int(ch) < 0x20
            @printf(out, "\\u%04x", Int(ch))
        else
            print(out, ch)
        end
    end
    String(take!(out))
end

function _json_value(x::AbstractString)::String
    "\"" * _json_escape(x) * "\""
end
_json_value(x::Bool)::String = x ? "true" : "false"
_json_value(x::Integer)::String = string(x)
_json_value(x::AbstractFloat)::String = isfinite(x) ? repr(x) : "null"
_json_value(::Nothing)::String = "null"

function _json_value(v::Vector)::String
    "[" * join([_json_value(x) for x in v], ", ") * "]"
end

function _json_value(d::Dict{String})::String
    pairs = [
        "\"" * _json_escape(k) * "\": " * _json_value(v) for (k, v) in sort(collect(d))
    ]
    "{" * join(pairs, ", ") * "}"
end

"""
    to_json_string(x) -> String

Serialise a nested structure of strings, numbers, booleans, `nothing`,
vectors and `Dict{String,...}` maps to JSON. Polynomials should be
pre-rendered with `lp_show_string`.
"""
to_json_string(x)::String = _json_value(x)

"""
    serialised_poly(p::LaurentPoly) -> String

Skein.jl-style serialisation of a Laurent polynomial: `exp:coeff` pairs
sorted by exponent and joined with commas.
"""
function serialised_poly(p::LaurentPoly)::String
    isempty(p) && return "0:0"
    join(["$(e):$(p[e])" for e in sort(collect(keys(p)))], ",")
end

"""
    export_knot_record(r::ConwaySuiteReport; gauss = Int[]) -> Dict{String}

A record shaped for Skein.jl's `KnotRecord` interchange: the same field
families (gauss code, writhe, serialised polynomials, determinant,
signature, genus, Seifert circle count) so downstream storage can ingest
KnotKnot evaluations without translation.
"""
function export_knot_record(r::ConwaySuiteReport; gauss::Vector{Int} = Int[])
    Dict{String,Any}(
        "name" => r.name,
        "gauss_code" => gauss,
        "diagram_format" => "knotknot-pd-v0",
        "crossing_number" => r.crossings,
        "writhe" => r.writhe,
        "alexander_polynomial" => serialised_poly(r.alexander),
        "jones_polynomial" => serialised_poly(r.jones),
        "conway_polynomial" => serialised_poly(r.conway),
        "determinant" => r.determinant,
        "signature" => r.signature,
        "genus" => r.genus,
        "seifert_circle_count" => r.seifert_circles,
        "colorings" => Dict{String,Any}(String(q) => c for (q, c) in r.colorings),
        "generator" => "KnotKnot.jl",
    )
end

"""
    write_visualization(path::AbstractString, entry::KnotEntry; extras = Dict{String,Any}())

Write the web-laboratory payload for one standard knot: invariants,
annotations and catalog metadata as JSON. The visualiser fetches this file
and layers theory annotations (quandle arcs, skein crossings, Conway core)
over the mechanics view.
"""
function write_visualization(
    path::AbstractString,
    entry::KnotEntry;
    extras::Dict{String,Any} = Dict{String,Any}(),
)::Nothing
    report = conway_suite(entry)
    payload = Dict{String,Any}(
        "schema" => "knotknot/visual-v0",
        "record" => export_knot_record(report; gauss = skein_gauss_code(entry.pd)),
        "pd" => [collect(t) for t in pd(entry.pd)],
        "description" => entry.description,
        "annotations" => Dict{String,Any}(
            "skein" => "Crossings carry the L+/L-/L0 skein triple; Conway's relation is verified per crossing.",
            "quandle" => "Arcs are fundamental-quandle generators; colouring counts are homomorphism counts.",
            "conway" => "Rational knots are classified by their Conway fraction p/q; determinant = |p|.",
        ),
        "extras" => extras,
    )
    open(path, "w") do io
        print(io, to_json_string(payload))
    end
    nothing
end

"""
    to_knottheory_dict(d::PlanarDiagram) -> Dict{String,Any}

KnotTheory.jl interop: the diagram as `entries` (KnotAtlas PD 5-tuples)
plus component grouping, ready for `KnotTheory.pdcode`.
"""
function to_knottheory_dict(d::PlanarDiagram)::Dict{String,Any}
    Dict{String,Any}(
        "entries" => [collect(t) for t in pd(d)],
        "components" => [collect(c) for c in d.components],
        "convention" => "knotatlas-pd: under a->c, over d->b, arcs counter-clockwise",
    )
end

"""
    from_knottheory_dict(x::Dict) -> PlanarDiagram

Build a `PlanarDiagram` from a KnotTheory.jl-style dict produced by
`to_knottheory_dict` (or by `KnotTheory.pdcode` round-tripping).
"""
function from_knottheory_dict(x::Dict)::PlanarDiagram
    entries = NTuple{5,Int}[]
    for row in x["entries"]
        push!(entries, (Int(row[1]), Int(row[2]), Int(row[3]), Int(row[4]), Int(row[5])))
    end
    comps = Vector{Int}[]
    if haskey(x, "components")
        for c in x["components"]
            push!(comps, Int[v for v in c])
        end
    end
    pdcode(entries; components = comps)
end

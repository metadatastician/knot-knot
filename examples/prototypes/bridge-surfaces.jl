# SPDX-License-Identifier: MPL-2.0
# bridge-surfaces.jl — minimal working prototype for the interop / report
# bridge surfaces (src/interop.jl + src/report.jl).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/bridge-surfaces.jl`; it exits
# non-zero unless every pinned assertion holds.
#
# What this proves: the three bridges to the outside world hold.
#
#   1. the JSON serialiser escapes what JSON requires, sorts object keys (so
#      output is byte-stable and hashable) and renders nothing / non-finite
#      floats as null;
#   2. the Skein.jl `exp:coeff` polynomial serialisation is exponent-sorted
#      and stable, with "0:0" for the zero polynomial;
#   3. the KnotTheory.jl PD dict round-trips exactly, crossings and component
#      grouping both, so a diagram can leave and come back unchanged.
#
# Conventions mirror src/interop.jl: PD 5-tuples in the KnotAtlas convention
# (under a->c, over d->b), the convention string
# "knotatlas-pd: under a->c, over d->b, arcs counter-clockwise", and the
# Skein.jl `KnotRecord` field families.

const FAILURES = Ref(0)

function check(label, got, want)
    if got == want
        println("  PASS  $label")
    else
        FAILURES[] += 1
        println("  FAIL  $label: got $got, want $want")
    end
end

# --- a minimal exact Laurent polynomial -------------------------------------

# Exponent -> coefficient, zero coefficients dropped so equality is structural.
const LaurentPoly = Dict{Int,Int}

lp_norm(p::LaurentPoly) = LaurentPoly(k => v for (k, v) in p if !iszero(v))

"""Skein.jl-style serialisation: `exp:coeff` pairs sorted by exponent."""
function serialised_poly(p::LaurentPoly)
    isempty(p) && return "0:0"
    join(["$(e):$(p[e])" for e in sort(collect(keys(p)))], ",")
end

# --- the JSON serialiser ----------------------------------------------------

function json_escape(s::AbstractString)
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
            print(out, "\\u" * lpad(string(Int(ch), base = 16), 4, '0'))
        else
            print(out, ch)
        end
    end
    String(take!(out))
end

json_value(x::AbstractString) = "\"" * json_escape(x) * "\""
json_value(x::Bool) = x ? "true" : "false"
json_value(x::Integer) = string(x)
json_value(x::AbstractFloat) = isfinite(x) ? repr(x) : "null"
json_value(::Nothing) = "null"

function json_value(v::Vector)
    "[" * join([json_value(x) for x in v], ", ") * "]"
end

function json_value(d::Dict{String})
    pairs = ["\"" * json_escape(k) * "\": " * json_value(v) for (k, v) in sort(collect(d))]
    "{" * join(pairs, ", ") * "}"
end

to_json_string(x) = json_value(x)

# --- a minimal planar diagram ----------------------------------------------

const Crossing = NTuple{5,Int}
const PlanarDiagram = Tuple{Vector{Crossing},Vector{Vector{Int}}}

"""PD 5-tuples in the KnotAtlas convention, plus component grouping."""
pd(d::PlanarDiagram) = d[1]

"""
    to_knottheory_dict(d)

KnotTheory.jl interop: the diagram as `entries` plus component grouping.
"""
function to_knottheory_dict(d::PlanarDiagram)
    Dict{String,Any}(
        "entries" => [collect(t) for t in pd(d)],
        "components" => [collect(c) for c in d[2]],
        "convention" => "knotatlas-pd: under a->c, over d->b, arcs counter-clockwise",
    )
end

"""Rebuild a diagram from a KnotTheory.jl-style dict."""
function from_knottheory_dict(x::Dict)
    entries = Crossing[]
    for row in x["entries"]
        push!(entries, (Int(row[1]), Int(row[2]), Int(row[3]), Int(row[4]), Int(row[5])))
    end
    comps = Vector{Int}[]
    if haskey(x, "components")
        for c in x["components"]
            push!(comps, Int[v for v in c])
        end
    end
    (entries, comps)
end

# --- the report / record bridge --------------------------------------------

"""The field families a ConwaySuiteReport carries (mirrored for the prototype)."""
struct SuiteReport
    name::String
    crossings::Int
    writhe::Int
    jones::LaurentPoly
    alexander::LaurentPoly
    conway::LaurentPoly
    determinant::Int
    signature::Int
    genus::Int
    seifert_circles::Int
    colorings::Vector{Pair{String,Int}}
end

"""A record shaped for Skein.jl's `KnotRecord` interchange."""
function export_knot_record(r::SuiteReport; gauss::Vector{Int} = Int[])
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

# --- the pins ---------------------------------------------------------------

const TREFOIL_PD = PlanarDiagram(
    [(1, 4, 2, 5, 1), (3, 6, 4, 1, 1), (5, 2, 6, 3, 1)],
    [[1, 2, 3, 4, 5, 6]],
)

const TREFOIL_REPORT = SuiteReport(
    "3_1", 3, 3,
    lp_norm(LaurentPoly(-4 => -1, -3 => 1, -1 => 1)),      # Jones
    lp_norm(LaurentPoly(-1 => 1, 0 => -1, 1 => 1)),        # Alexander
    lp_norm(LaurentPoly(0 => 1, 2 => 1)),                  # Conway
    3, -2, 1, 2,
    ["R3" => 9, "R5" => 5],
)

function main()
    println("bridge-surfaces prototype")
    println("  (JSON serialiser, Skein.jl polynomial records, KnotTheory PD round-trip)")

    println("\n[1] polynomial serialisation: exp:coeff, exponent-sorted")
    check("trefoil Alexander", serialised_poly(lp_norm(LaurentPoly(-1 => 1, 0 => -1, 1 => 1))),
          "-1:1,0:-1,1:1")
    check("the zero polynomial is \"0:0\"", serialised_poly(lp_norm(LaurentPoly())), "0:0")
    check("a single constant term", serialised_poly(lp_norm(LaurentPoly(0 => 5))), "0:5")
    check("negative coefficients survive", serialised_poly(lp_norm(LaurentPoly(0 => 1, 2 => -1))),
          "0:1,2:-1")
    check("gaps in the exponents are preserved",
          serialised_poly(lp_norm(LaurentPoly(-3 => 1, 0 => -1, 7 => 2))), "-3:1,0:-1,7:2")
    check("trefoil Conway", serialised_poly(lp_norm(LaurentPoly(0 => 1, 2 => 1))), "0:1,2:1")
    check("figure-eight Alexander",
          serialised_poly(lp_norm(LaurentPoly(-1 => -1, 0 => 3, 1 => -1))), "-1:-1,0:3,1:-1")
    check("cinquefoil Alexander",
          serialised_poly(lp_norm(LaurentPoly(-2 => 1, -1 => -1, 0 => 1, 1 => -1, 2 => 1))),
          "-2:1,-1:-1,0:1,1:-1,2:1")

    println("\n[2] JSON escaping")
    check("a double quote", json_escape("a\"b"), "a\\\"b")
    check("a backslash", json_escape("a\\b"), "a\\\\b")
    check("a newline", json_escape("a\nb"), "a\\nb")
    check("a tab", json_escape("a\tb"), "a\\tb")
    check("a carriage return", json_escape("a\rb"), "a\\rb")
    check("a control character", json_escape("a\x01b"), "a\\u0001b")
    check("the highest control character", json_escape("\x1f"), "\\u001f")
    check("a space is not escaped", json_escape("a b"), "a b")
    check("non-ASCII passes through", json_escape("éè"), "éè")
    check("an empty string", json_escape(""), "")

    println("\n[3] JSON value rendering")
    check("strings are quoted", json_value("x"), "\"x\"")
    check("true and false", json_value(true), "true")
    check("nothing is null", json_value(nothing), "null")
    check("a non-finite float is null", json_value(Inf), "null")
    check("NaN is null", json_value(NaN), "null")
    check("a finite float", json_value(1.5), "1.5")
    check("integers", json_value(42), "42")
    check("a nested vector", json_value([1, "a", nothing]), "[1, \"a\", null]")
    check("an empty vector", json_value(Int[]), "[]")
    check("an empty object", json_value(Dict{String,Any}()), "{}")
    # keys are sorted, so the output is byte-stable
    check("object keys are sorted",
          json_value(Dict{String,Any}("b" => 1, "a" => 2, "c" => 3)),
          "{\"a\": 2, \"b\": 1, \"c\": 3}")
    check("nested objects sort at every level",
          json_value(Dict{String,Any}("z" => Dict{String,Any}("y" => 1, "x" => 2))),
          "{\"z\": {\"x\": 2, \"y\": 1}}")

    println("\n[4] the KnotTheory PD round-trip")
    dict = to_knottheory_dict(TREFOIL_PD)
    check("entries are the PD 5-tuples",
          dict["entries"], [[1, 4, 2, 5, 1], [3, 6, 4, 1, 1], [5, 2, 6, 3, 1]])
    check("components are preserved", dict["components"], [[1, 2, 3, 4, 5, 6]])
    check("the convention string is stable",
          dict["convention"], "knotatlas-pd: under a->c, over d->b, arcs counter-clockwise")
    back = from_knottheory_dict(dict)
    check("crossings survive the round trip", pd(back), pd(TREFOIL_PD))
    check("components survive the round trip", back[2], TREFOIL_PD[2])
    check("the round trip is an involution",
          to_knottheory_dict(from_knottheory_dict(to_knottheory_dict(TREFOIL_PD))),
          to_knottheory_dict(TREFOIL_PD))
    # a multi-component diagram keeps its grouping
    hopf = PlanarDiagram([(1, 3, 2, 4, 1), (2, 4, 1, 3, 1)], [[1, 2], [3, 4]])
    hdict = to_knottheory_dict(hopf)
    check("a two-component diagram keeps both components",
          from_knottheory_dict(hdict)[2], [[1, 2], [3, 4]])
    # a missing components key must not lose the crossings
    stripped = Dict{String,Any}("entries" => dict["entries"])
    check("a dict without components still rebuilds the crossings",
          pd(from_knottheory_dict(stripped)), pd(TREFOIL_PD))
    check("and reports no components", from_knottheory_dict(stripped)[2], Vector{Int}[])

    println("\n[5] the Skein.jl record")
    rec = export_knot_record(TREFOIL_REPORT; gauss = [-1, 3, -2, 1, -3, 2])
    check("the field families are exactly these",
          sort(collect(keys(rec))),
          sort(["alexander_polynomial", "colorings", "crossing_number", "determinant",
                "diagram_format", "gauss_code", "generator", "genus", "jones_polynomial",
                "name", "conway_polynomial", "seifert_circle_count", "signature", "writhe"]))
    check("the diagram format tag", rec["diagram_format"], "knotknot-pd-v0")
    check("the generator tag", rec["generator"], "KnotKnot.jl")
    check("the gauss code is passed through", rec["gauss_code"], [-1, 3, -2, 1, -3, 2])
    check("the serialised Alexander", rec["alexander_polynomial"], "-1:1,0:-1,1:1")
    check("the serialised Jones", rec["jones_polynomial"], "-4:-1,-3:1,-1:1")
    check("the serialised Conway", rec["conway_polynomial"], "0:1,2:1")
    check("the invariants", (rec["determinant"], rec["signature"], rec["genus"]), (3, -2, 1))
    check("the colouring counts", rec["colorings"], Dict{String,Any}("R3" => 9, "R5" => 5))
    check("the default gauss code is empty",
          export_knot_record(TREFOIL_REPORT)["gauss_code"], Int[])
    # the record must actually be serialisable
    text = to_json_string(rec)
    check("the record serialises to JSON without error", occursin("\"determinant\": 3", text), true)
    check("the serialised record is sorted and stable",
          to_json_string(export_knot_record(TREFOIL_REPORT; gauss = Int[])),
          to_json_string(export_knot_record(TREFOIL_REPORT; gauss = Int[])))

    println("\n[6] structural sanity of the emitted JSON")
    # A light structural validator: balanced delimiters, no trailing commas,
    # every key quoted. Not a parser, but it catches the usual mistakes.
    function json_well_formed(s::AbstractString)
        depth = 0
        in_string = false
        escaped = false
        for ch in s
            if escaped
                escaped = false
                continue
            end
            if ch == '\\'
                escaped = in_string
                continue
            end
            if ch == '"'
                in_string = !in_string
                continue
            end
            in_string && continue
            if ch in ('{', '[')
                depth += 1
            elseif ch in ('}', ']')
                depth -= 1
                depth < 0 && return false
            end
        end
        (depth == 0 && !in_string) || return false
        !occursin(r",\s*[}\]]", s) || return false
        !occursin(r"\{\s*,", s) || return false
        true
    end
    for (label, payload) in (
            ("the knot record", rec),
            ("the PD dict", dict),
            ("a nested payload",
             Dict{String,Any}("record" => rec, "pd" => [collect(t) for t in pd(TREFOIL_PD)],
                              "extras" => Dict{String,Any}("note" => "a \"quoted\" note"))))
        check("$label is well-formed JSON", json_well_formed(to_json_string(payload)), true)
    end

    println()
    if FAILURES[] == 0
        println("ALL PASS — the interop/report bridge-surface prototype is sound.")
    else
        println("$(FAILURES[]) FAILURE(S) — the interop/report prototype is wrong.")
        exit(1)
    end
end

main()

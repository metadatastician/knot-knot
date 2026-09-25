# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The Conway-suite report: one call evaluates everything the suite needs
# for a diagram and renders it for humans and machines.

"""
    ConwaySuiteReport

Full invariant suite for one knot diagram: crossing data, Jones /
Alexander / Conway polynomials, determinant, signature, genus, Seifert
circle count, and finite-quandle colouring counts.
"""
struct ConwaySuiteReport
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

"""
    conway_suite(x; name = "") -> ConwaySuiteReport

Evaluate the Conway suite on a `PlanarDiagram`, a `KnotEntry`, or a
standard-knot name (`"3_1"`...). Colouring counts are evaluated against
the dihedral quandles R3 and R5.
"""
function conway_suite(d::PlanarDiagram; name::String = "")::ConwaySuiteReport
    colorings = Pair{String,Int}[
        "R3" => coloring_count(d, dihedral_quandle(3)),
        "R5" => coloring_count(d, dihedral_quandle(5)),
    ]
    ConwaySuiteReport(
        name,
        crossing_count(d),
        writhe(d),
        jones_polynomial(d),
        alexander_polynomial(d),
        conway_polynomial(d),
        determinant(d),
        signature(d),
        genus(d),
        seifert_circles(d),
        colorings,
    )
end

function conway_suite(entry::KnotEntry; name::String = entry.name)::ConwaySuiteReport
    conway_suite(entry.pd; name = name)
end

function conway_suite(name::AbstractString)::ConwaySuiteReport
    conway_suite(standard_knot(name))
end

function Base.show(io::IO, r::ConwaySuiteReport)
    label = isempty(r.name) ? "diagram" : r.name
    println(io, "Conway suite — $(label)")
    println(io, "  crossings        : $(r.crossings)   writhe: $(r.writhe)")
    println(io, "  Jones V(t)       : $(lp_show_string(r.jones))")
    println(io, "  Alexander D(t)   : $(lp_show_string(r.alexander))")
    println(io, "  Conway nabla(z)  : $(lp_show_string(r.conway; var = "z"))")
    println(io, "  determinant      : $(r.determinant)   signature: $(r.signature)")
    println(io, "  genus            : $(r.genus)   Seifert circles: $(r.seifert_circles)")
    for (q, c) in r.colorings
        println(io, "  $(q) colourings    : $(c)")
    end
end

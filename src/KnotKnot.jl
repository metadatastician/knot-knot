# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# KnotKnot.jl — knot theory across the theory-practice divide.
#
# Theory side: diagram kernels (planar / Gauss codes), Jones and
# Alexander-Conway invariants, quandles, braid closures, rational tangle
# (Conway) calculus, and skein-relation machinery.
#
# Practice side: elastic-rod knot mechanics (capstan friction, bend-stress
# superposition, rope-efficiency prediction) and protistological topology
# (flagellar waveforms, writhe/linking, kinetoplast catenation).
#
# Bridge surfaces: the Conway-suite report, JSON export for the web
# laboratory, and interop with the estate engines KnotTheory.jl / Skein.jl.
#
# Minted from the hyperpolymath julia-library archetype; provenance in
# .machine_readable/PROVENANCE.a2ml. Conventions are pinned against
# KnotTheory.jl (KnotAtlas PD codes) and Skein.jl (signed Gauss codes);
# see docs/theory/mathematics/KNOT-THEORY.adoc.

module KnotKnot

using LinearAlgebra
using Printf

# --- Laurent polynomials (integer coefficients) -------------------------
include("laurent.jl")

# --- Diagram kernel: PD codes, Gauss traversals, standard table ---------
include("diagrams.jl")
include("gauss.jl")
include("knot_table.jl")

# --- Invariants ----------------------------------------------------------
include("alexander.jl")
include("jones.jl")
include("seifert.jl")

# --- Algebraic structures: quandles, braids, Conway tangles, skein ------
include("quandle.jl")
include("braid.jl")
include("tangle.jl")
include("skein.jl")

# --- Applied layers -------------------------------------------------------
include("mechanics.jl")
include("protistology.jl")

# --- Bridge surfaces ------------------------------------------------------
include("report.jl")
include("interop.jl")

# --- public API -----------------------------------------------------------
export LaurentPoly, lp_eval, lp_show_string, lp_pow, lp_shift
export Crossing, PlanarDiagram, pdcode, pd, writhe, mirror, crossing_count
export arcs_of, traverse_components
export GaussWord, gauss_word, skein_gauss_code
export standard_knot, knot_names, KnotEntry
export alexander_polynomial, conway_polynomial, determinant, genus
export jones_polynomial, bracket_polynomial
export seifert_circles, seifert_matrix, signature
export Quandle, dihedral_quandle, alexander_quandle, is_quandle
export fundamental_relations, coloring_count
export Braid, braid_closure, parse_braid_word
export RationalTangle, tangle_fraction, conway_notation, continued_fraction, rotate
export two_bridge_det, is_knot_fraction
export switch_crossing, smooth_crossing, skein_triple, verify_conway_skein
export Rod, curvature_profile, bending_energy, capstan_tensions
export stress_profile, safety_factor, breaking_force, knot_efficiency
export writhe_of, linking_of, calugareanu_check
export flagellum, beat_summary, kinetoplast_chain, protistology_report, linking_matrix
export ConwaySuiteReport, conway_suite
export to_json_string, export_knot_record, write_visualization
export to_knottheory_dict, from_knottheory_dict

end # module KnotKnot

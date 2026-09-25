# SPDX-License-Identifier: MPL-2.0
# make.jl — Documenter build for KnotKnot.jl.
#
# Deploy posture: the estate workflow builds these pages on every push to
# main; deployment requires a DOCUMENTER_KEY deploy key that the repository
# does not hold yet (docs/decisions/0002-docs-build-only.adoc). The build
# step is wired and green either way.

using Documenter
using KnotKnot

makedocs(
    sitename = "KnotKnot.jl",
    modules = [KnotKnot],
    pages = [
        "Home" => "index.adoc",
        "Theory" => "theory.adoc",
        "Conway suite" => "conway.adoc",
        "Mechanics" => "mechanics.adoc",
        "Protistology" => "protistology.adoc",
    ],
    warnonly = true,
)

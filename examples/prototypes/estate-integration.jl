# SPDX-License-Identifier: MPL-2.0
# estate-integration.jl — minimal working prototype for ULTRAPLAN Phase 9
# (estate integration).
#
# Self-contained: Julia standard library only, no KnotKnot import, no
# Project.toml. Run `julia examples/prototypes/estate-integration.jl`; it exits
# non-zero unless every pinned assertion holds.
#
# Phase 9 is not a mathematics kernel, so the prototype encodes the RULES the
# integration must satisfy rather than a computation:
#
#   1. the gv-clade-index CLADE record is internally consistent — the closed
#      taxonomy of twelve clades, the code/name agreement CLADE-006 checks,
#      the prefixed name, the uuid format, the lineage and status taxonomies;
#   2. the General-registry metadata is well-formed — package name, uuid,
#      semver, git-tree-sha1 shape and the Julia compat bound;
#   3. the estate lockfile rule holds — a STEP-level `uses:` must be recorded
#      in actions.lock under that workflow's own path, while a JOB-level
#      reusable-workflow miss is harmless. This is the rule whose violation
#      leaves 39 estate repos with silently dead CI (jobs = 0,
#      startup_failure).
#
# The values below mirror .machine_readable/descriptiles/CLADE.a2ml and
# Project.toml as they stand; they are embedded so the prototype runs anywhere.

const FAILURES = Ref(0)

function check(label, got, want)
    if got == want
        println("  PASS  $label")
    else
        FAILURES[] += 1
        println("  FAIL  $label: got $got, want $want")
    end
end

# --- the closed clade taxonomy ----------------------------------------------

const CLADE_NAMES = Dict{String,String}(
    "fv" => "Formal Verification & Proofs",
    "nl" => "Nextgen Languages",
    "rm" => "Repo Management & Tooling",
    "gv" => "Governance & Standards",
    "db" => "Databases",
    "ap" => "Applications",
    "ix" => "Infrastructure & Cloud",
    "dx" => "Developer Ecosystem",
    "pt" => "Protocols & Interop",
    "ax" => "AI & Neurosymbolic",
    "gm" => "Games & Interactive",
    "sc" => "Security",
)

const LINEAGE_TYPES = Set([
    "standalone", "monorepo", "monorepo-child", "inflated", "deflated", "hub", "satellite",
])

const STATUS_PHASES = Set([
    "reserved", "incubating", "active", "dormant", "merged", "superseded", "archived", "extinct",
])

# Phases for which the repo is expected to still exist.
const LIVE_PHASES = Set(["reserved", "incubating", "active", "dormant"])

# --- the record as it stands ------------------------------------------------

const CLADE = Dict{String,Any}(
    "uuid" => "a6cbd9d9-f370-519d-a0d0-9d00bf3e780b",
    "primary-forge" => "github",
    "primary-owner" => "metadatastician",
    "canonical-name" => "knot-knot",
    "prefixed-name" => "dx-knot-knot",
    "clade" => Dict{String,Any}(
        "primary" => "dx",
        "primary-name" => "Developer Ecosystem",
        "secondary" => ["fv"],
        "assigned" => "2026-09-25",
        "rationale" => "Value proposition is a Julia toolkit other science and engineering work builds on.",
    ),
    "forges" => Dict{String,Any}("github" => "metadatastician/knot-knot",
                                 "gitlab" => "", "bitbucket" => ""),
    "lineage" => Dict{String,Any}("type" => "standalone", "parent" => "", "born" => "2026-09-25",
                                  "instantiated-from" => "rsr-template-repo"),
    "status" => Dict{String,Any}("phase" => "incubating", "since" => "2026-09-25", "present" => true),
)

const PROJECT = Dict{String,Any}(
    "name" => "KnotKnot",
    "uuid" => "74a22736-c88b-55db-9726-7f94661135a1",
    "version" => "0.1.0",
    "repo" => "https://github.com/metadatastician/knot-knot.git",
    "compat_julia" => "1.9",
    "deps" => Dict{String,String}(
        "LinearAlgebra" => "37e2e46d-f89d-539d-b4ee-838fcccc9c8e",
        "Printf" => "de0858da-6303-5e67-8744-51eddeeeb8d7",
    ),
    "extras" => Dict{String,String}(
        "Test" => "8dfed614-e22c-5e08-85e1-65c5234f0b40",
        "Aqua" => "4c88cf16-eb10-579e-8560-4a9242c79595",
    ),
)

# --- validation helpers -----------------------------------------------------

function is_uuid(s::AbstractString)
    m = match(r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", s)
    m === nothing && return false
    # the variant nibble (first character of the fourth group) must be 8, 9, a or b
    parts = split(s, "-")
    parts[4][1] in ('8', '9', 'a', 'b') || return false
    true
end

"""UUID version nibble: the first character of the third group."""
uuid_version(s::AbstractString) = parse(Int, split(s, "-")[3][1:1], base = 16)

"""Semver `major.minor.patch`, no prerelease or build metadata."""
is_semver(s::AbstractString) =
    match(r"^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$", s) !== nothing

"""A git tree object id: 40 lowercase hex characters."""
is_tree_sha1(s::AbstractString) =
    match(r"^[0-9a-f]{40}$", s) !== nothing

# --- the pins ---------------------------------------------------------------

function main()
    println("estate-integration prototype")
    println("  (clade record, registry readiness, estate lockfile rule)")

    println("\n[1] the clade taxonomy is closed and correctly named")
    check("there are exactly twelve clades", length(CLADE_NAMES), 12)
    check("every code is two letters",
          all(length(c) == 2 for c in keys(CLADE_NAMES)), true)
    check("dx is Developer Ecosystem", CLADE_NAMES["dx"], "Developer Ecosystem")
    check("fv is Formal Verification & Proofs", CLADE_NAMES["fv"], "Formal Verification & Proofs")
    check("pt is Protocols & Interop (the paint-type trap)",
          CLADE_NAMES["pt"], "Protocols & Interop")
    check("gv is Governance & Standards (the gossamer trap)",
          CLADE_NAMES["gv"], "Governance & Standards")
    check("rm is Repo Management & Tooling (the template default)",
          CLADE_NAMES["rm"], "Repo Management & Tooling")

    println("\n[2] this repo's clade record is internally consistent")
    check("primary is in the taxonomy", haskey(CLADE_NAMES, CLADE["clade"]["primary"]), true)
    check("primary-name agrees with primary (CLADE-006)",
          CLADE["clade"]["primary-name"], CLADE_NAMES[CLADE["clade"]["primary"]])
    check("prefixed-name is primary + '-' + canonical-name",
          CLADE["prefixed-name"], CLADE["clade"]["primary"] * "-" * CLADE["canonical-name"])
    check("secondary codes are in the taxonomy",
          all(haskey(CLADE_NAMES, c) for c in CLADE["clade"]["secondary"]), true)
    check("secondary does not repeat primary",
          !(CLADE["clade"]["primary"] in CLADE["clade"]["secondary"]), true)
    check("the clade uuid is well-formed", is_uuid(CLADE["uuid"]), true)
    check("the clade uuid is version 5 (derived, not random)",
          uuid_version(CLADE["uuid"]), 5)
    check("the primary forge is github", CLADE["primary-forge"], "github")
    check("the primary owner is metadatastician", CLADE["primary-owner"], "metadatastician")
    check("the github forge path is owner/name",
          CLADE["forges"]["github"], "metadatastician/knot-knot")
    check("the unpopulated forges are empty strings",
          (CLADE["forges"]["gitlab"], CLADE["forges"]["bitbucket"]), ("", ""))

    println("\n[3] lineage and status taxonomies")
    check("the lineage type is in the closed set",
          CLADE["lineage"]["type"] in LINEAGE_TYPES, true)
    check("the status phase is in the closed set",
          CLADE["status"]["phase"] in STATUS_PHASES, true)
    check("a live phase has present = true",
          CLADE["status"]["present"], CLADE["status"]["phase"] in LIVE_PHASES)
    check("the assignment date is ISO", match(r"^\d{4}-\d{2}-\d{2}$", CLADE["clade"]["assigned"]) !== nothing, true)
    check("the birth date is ISO", match(r"^\d{4}-\d{2}-\d{2}$", CLADE["lineage"]["born"]) !== nothing, true)

    println("\n[4] General-registry readiness")
    check("the package name matches the module name", PROJECT["name"], "KnotKnot")
    check("the registry uuid is well-formed", is_uuid(PROJECT["uuid"]), true)
    check("the package uuid differs from the clade uuid",
          PROJECT["uuid"] != CLADE["uuid"], true)
    check("the version is semver", is_semver(PROJECT["version"]), true)
    check("the compat bound is a bare minor version",
          match(r"^[0-9]+\.[0-9]+$", PROJECT["compat_julia"]) !== nothing, true)
    check("the repo URL is the canonical forge URL",
          PROJECT["repo"], "https://github.com/metadatastician/knot-knot.git")
    check("every dependency has a uuid",
          all(is_uuid(v) for v in values(PROJECT["deps"])), true)
    check("every extra has a uuid",
          all(is_uuid(v) for v in values(PROJECT["extras"])), true)
    check("the two stdlib deps are LinearAlgebra and Printf",
          sort(collect(keys(PROJECT["deps"]))), ["LinearAlgebra", "Printf"])
    # a tree-sha1 is what the registry records per version
    check("a git tree object id is 40 lowercase hex characters",
          is_tree_sha1("0123456789abcdef0123456789abcdef01234567"), true)
    check("an uppercase tree object id is rejected",
          is_tree_sha1("0123456789ABCDEF0123456789ABCDEF01234567"), false)
    check("a short tree object id is rejected", is_tree_sha1("0123456789abcdef"), false)
    check("a non-hex tree object id is rejected",
          is_tree_sha1("z123456789abcdef0123456789abcdef01234567"), false)

    println("\n[5] the estate lockfile rule")
    # Rule (measured on 39 estate repos): a STEP-level `uses:` that is absent
    # from actions.lock under that workflow's own path kills the run at
    # startup (jobs = 0, startup_failure). A JOB-level reusable-workflow ref
    # that is missing is harmless. Matching strips the subpath and is
    # case-insensitive.
    LOCK = Dict{String,Any}(
        ".github/workflows/ci.yml" => [
            "actions/checkout@v4", "julia-actions/setup-julia@v2",
            "julia-actions/cache@v2", "julia-actions/julia-buildpkg@v1",
            "julia-actions/julia-runtest@v1",
        ],
        ".github/workflows/docs.yml" => [
            "actions/checkout@v4", "julia-actions/setup-julia@v2",
            "julia-actions/julia-docdeploy@v1",
        ],
        ".github/workflows/release.yml" => [
            "actions/checkout@v4", "julia-actions/setup-julia@v2",
        ],
    )

    """Strip the subpath and lowercase: `owner/repo/sub/path@ref` -> `owner/repo@ref`."""
    function lock_key(ref::AbstractString)
        at = findlast('@', ref)
        at === nothing && return lowercase(ref)
        i = last(at)
        path = ref[1:(i - 1)]
        parts = split(path, '/')
        owner_repo = length(parts) >= 2 ? parts[1] * "/" * parts[2] : path
        lowercase(owner_repo * ref[i:end])
    end

    function step_recorded(workflow::AbstractString, ref::AbstractString)
        entries = get(LOCK, workflow, String[])
        any(e -> lock_key(e) == lock_key(ref), entries)
    end

    # every step-level ref in the pinned workflows is recorded
    for (wf, refs) in sort(collect(LOCK))
        check("$wf: every step-level uses is locked under its own path",
              all(step_recorded(wf, r) for r in refs), true)
    end
    # a subpath ref still matches its owner/repo entry
    check("a subpath ref matches the owner/repo lock entry",
          step_recorded(".github/workflows/ci.yml", "actions/checkout/sub/path@v4"), true)
    check("matching is case-insensitive",
          step_recorded(".github/workflows/ci.yml", "Actions/CheckOut@V4"), true)
    # a ref locked under a DIFFERENT workflow's path does not count
    check("a ref locked only under another workflow is not resolvable",
          step_recorded(".github/workflows/ci.yml", "julia-actions/julia-docdeploy@v1"), false)
    # an unknown ref is not resolvable
    check("an unknown ref is not resolvable",
          step_recorded(".github/workflows/ci.yml", "actions/upload-artifact@v4"), false)
    # a workflow with no lock entries at all resolves nothing
    check("a workflow with no lock entries resolves nothing",
          step_recorded(".github/workflows/new.yml", "actions/checkout@v4"), false)

    println("\n[6] what Phase 9 still blocks on")
    # Registration is the owner's act: deriving the uuid does not register the
    # repo, and an agent must not do it.
    check("the clade is chosen (owner ratification pending)", CLADE["clade"]["primary"], "dx")
    check("the repo is not yet registered: the status is incubating",
          CLADE["status"]["phase"], "incubating")
    check("derivation is not registration: the phase is still incubating",
          CLADE["status"]["phase"], "incubating")
    # the estate engines this package bridges to
    ENGINES = ["KnotTheory.jl", "Skein.jl", "Axiom.jl"]
    check("the bridged estate engines", sort(ENGINES), ["Axiom.jl", "KnotTheory.jl", "Skein.jl"])

    println()
    if FAILURES[] == 0
        println("ALL PASS — the estate-integration prototype is sound.")
    else
        println("$(FAILURES[]) FAILURE(S) — the estate-integration prototype is wrong.")
        exit(1)
    end
end

main()

# SPDX-License-Identifier: MPL-2.0
# axiom_bridge.jl — hand-off surface from KnotKnot mechanics to Axiom.jl.
#
# Axiom.jl (the estate's ML / proof framework) learns from exported data;
# it is deliberately NOT a dependency of KnotKnot. This example builds the
# training set Axiom would consume: (curvature, local tension) -> stress,
# sampled along a knotted rod, and writes it as JSON plus a small proof
# obligation stub (the mechanics identities Axiom should witness).

using KnotKnot

function knotted_rod(; samples = 240)
    pts = Vector{NTuple{3,Float64}}()
    for i in 0:samples
        t = i / samples
        x = -14.0 + 28.0 * t
        bump = exp(-((t - 0.5)^2) / 0.005)
        y = bump * 2.2 * sin(6.0 * pi * t)
        z = bump * 2.2 * cos(6.0 * pi * t) * 0.4
        push!(pts, (x, y, z))
    end
    Rod(pts; diameter = 2.0, youngs_mpa = 3500.0, uts_mpa = 850.0, friction = 0.35)
end

rod = knotted_rod()
kappa = curvature_profile(rod)
sigma = stress_profile(rod, 180.0)

pairs = [
    Dict{String,Any}(
        "curvature" => kappa[i],
        "tension" => 180.0,
        "sigma_axial" => sigma.axial[i],
        "sigma_bending" => sigma.bending[i],
        "sigma_total" => sigma.total[i],
    ) for i in 1:2:length(kappa)
]

payload = Dict{String,Any}(
    "schema" => "knotknot/axiom-training-v0",
    "target" => "Axiom.jl regression: curvature x tension -> sigma_total",
    "features" => ["curvature", "tension"],
    "labels" => ["sigma_total"],
    "rows" => pairs,
    "proof_obligations" => [
        "sigma_total == sigma_axial + sigma_bending (superposition)",
        "sigma_bending == E * kappa * d / 2 (bend-stress identity)",
        "F(s) == F0 * exp(-mu * theta(s)) (capstan law)",
    ],
    "invariants" => Dict{String,Any}(
        "safety_factor" => safety_factor(rod, 180.0),
        "predicted_breaking_force" => breaking_force(rod),
        "efficiency" => knot_efficiency(rod),
    ),
)

path = joinpath(@__DIR__, "axiom_training_data.json")
open(path, "w") do io
    print(io, to_json_string(payload))
end
println("wrote $(path): $(length(pairs)) training rows + 3 proof obligations")

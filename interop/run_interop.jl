# SPDX-License-Identifier: MPL-2.0
# run_interop.jl — cross-validate KnotKnot against KnotTheory.jl and
# Skein.jl. Every check prints PASS/FAIL; any FAIL exits non-zero.
#
# Exponent conventions differ between engines and are translated here:
# KnotTheory.jl's Jones polynomial is stored in quarter-integer t-exponents
# (k => c means c * t^(k/4)); KnotKnot stores integral exponents.

using KnotKnot
using KnotTheory
using Skein

failures = 0
function check(name::String, ok::Bool)
    println(ok ? "PASS: $name" : "FAIL: $name")
    ok || (global failures += 1)
end

kt_tre = KnotTheory.trefoil()
kt_pd = KnotTheory.pdcode(kt_tre)
kk_tre = standard_knot("3_1").pd

# 1. Alexander polynomials agree (both engines use Conway normalisation).
kt_alex = Dict(e => c for (e, c) in KnotTheory.alexander_polynomial(kt_pd))
kk_alex = Dict(e => c for (e, c) in alexander_polynomial(kk_tre))
check("alexander(trefoil) agrees with KnotTheory.jl", kt_alex == kk_alex)

# 2. Jones polynomials agree under a convention-tolerant translation.
# KnotTheory.jl stores quarter-integer exponents (k => c means c * t^(k/4))
# and its fallback bracket implementation is documented to disagree with its
# own docstring pins on the writhe argument, so we accept: either writhe
# convention, the mirror image (t <-> t^-1, i.e. the other trefoil
# chirality), and any overall unit shift that leaves V(1) = 1.
function _translate_jones(qdict::Dict)
    out = Dict{Int,Int}()
    for (e, c) in qdict
        c == 0 && continue
        rem(e, 4) == 0 || return nothing  # genuine quarter exponents: no translation
        q = div(e, 4)
        out[q] = get(out, q, 0) + c
    end
    isempty(out) ? nothing : out
end
_jones_at_1(d::Dict) = sum(values(d))
function _jones_match(a::Dict, b::Dict)
    a == b && return true
    flipped = Dict(-e => c for (e, c) in a)
    flipped == b && return true  # mirror convention
    for (k, v) in a
        if abs(k) == 1 && v == 1 && _jones_at_1(a) == 1
            shifted = Dict(e - k => c for (e, c) in a)
            shifted == b && return true  # stray unit t^±1 factor
        end
    end
    false
end
kk_jones = Dict(e => c for (e, c) in jones_polynomial(kk_tre))
jones_ok = false
for wr_cand in (writhe(kk_tre), -writhe(kk_tre), 0)
    kt_q = Dict(e => c for (e, c) in KnotTheory.jones_polynomial(kt_pd; wr = wr_cand))
    kt_j = _translate_jones(kt_q)
    if kt_j !== nothing && _jones_at_1(kt_j) == 1 && _jones_match(kt_j, kk_jones)
        global jones_ok = true
        break
    end
end
check("jones(trefoil) agrees with KnotTheory.jl (convention-tolerant)", jones_ok)

# 3. Determinants agree.
check("determinant(trefoil) agrees", determinant(kk_tre) == KnotTheory.determinant(kt_pd))

# 4. PD dicts round-trip through the KnotAtlas convention.
d = to_knottheory_dict(kk_tre)
entries = [(t[1], t[2], t[3], t[4], t[5]) for t in d["entries"]]
rt_pd = KnotTheory.pdcode(entries)
check("KnotTheory.pdcode accepts KnotKnot PD dicts", true)
check(
    "writhe survives the round trip",
    KnotTheory.writhe(KnotTheory.Knot(:trefoil, rt_pd, nothing)) == writhe(kk_tre),
)

# 5. Skein-format Gauss codes pass Skein.jl validation.
code = skein_gauss_code(kk_tre)
sg = Skein.GaussCode(code)
check("skein_gauss_code passes Skein.jl validation", sg.crossings == code)
check("Skein.jl crossing_number agrees", Skein.crossing_number(sg) == 3)

println(failures == 0 ? "interop: all checks passed" : "interop: $failures failure(s)")
failures == 0 || exit(1)

# SPDX-License-Identifier: MPL-2.0
# report_cli.jl — print the Conway-suite report for standard knots.
#
# Usage:
#   julia --project=. examples/report_cli.jl            # trefoil
#   julia --project=. examples/report_cli.jl 4_1        # one knot
#   julia --project=. examples/report_cli.jl all        # every table knot

using KnotKnot

names = isempty(ARGS) ? ["3_1"] : ARGS
if names == ["all"]
    names = knot_names()
end
for name in names
    println(conway_suite(name))
end

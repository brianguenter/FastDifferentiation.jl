#!/usr/bin/env julia
# Checks `jacobian` on random expression graphs of several sizes and reports
# the time per graph. Any failure is printed as Julia code that reproduces it.
#
# Usage (from the repository root):
#     julia --project=. scripts/fuzz_factoring.jl [seed] [seconds_per_size]
#
# Defaults: seed 1, 20 seconds per size. Each size runs graphs until its time
# budget is used up, so total run time is about (number of sizes) x budget,
# plus compilation.

include(joinpath(@__DIR__, "..", "test", "random_graphs.jl"))
using .RandomGraphs
import FastDifferentiation as FD
using Random

seed = parse(Int, get(ARGS, 1, "1"))
budget = parse(Float64, get(ARGS, 2, "20"))

# (label, nvars range, nops range, nroots range)
const SIZES = [
    ("tiny",    2:3,   3:6,     1:3),
    ("small",   2:5,   8:15,    1:4),
    ("medium",  3:8,   20:40,   2:6),
    ("large",   4:12,  60:120,  3:10),
    ("x-large", 8:20,  200:400, 5:20),
]

rng = MersenneTwister(seed)
RandomGraphs.fuzz(rng, 20; nvars=2:3, nops=3:6, nroots=1:3)   # compile before timing

println(rpad("size", 9), rpad("vars", 7), rpad("ops", 9), rpad("roots", 7), lpad("graphs", 8), lpad("ms/graph", 10), lpad("skipped", 9), lpad("failures", 10))
total_failures = 0
for (label, nv, no, nr) in SIZES
    n = 0; skipped = 0; failures = Any[]
    t0 = time()
    while time() - t0 < budget
        f, s = RandomGraphs.fuzz(rng, 1; nvars=nv, nops=no, nroots=nr)
        append!(failures, f); skipped += s
        n += 1
        n % 200 == 0 && empty!(FD.EXPRESSION_CACHE)  # keep the global expression cache from growing without bound
    end
    elapsed = time() - t0
    println(rpad(label, 9), rpad("$(first(nv))-$(last(nv))", 7), rpad("$(first(no))-$(last(no))", 9), rpad("$(first(nr))-$(last(nr))", 7),
        lpad(n, 8), lpad(round(1000 * elapsed / n, digits=2), 10), lpad(skipped, 9), lpad(length(failures), 10))
    for (spec, status, msg) in first(failures, 3)
        println("\n  $status:\n", replace(spec_code(spec), r"^"m => "    "), "\n  ", first(msg, 400), "\n")
    end
    global total_failures += length(failures)
end
exit(total_failures == 0 ? 0 : 1)

# Summary statistics quoted in the report, computed from data/results.json.
#   julia investigations/tracingjit/report/stats.jl
using Statistics
txt = read(joinpath(@__DIR__, "data", "results.json"), String)
keys_ = ["status", "name", "pkg", "err", "speed", "hit", "fd", "first", "guards", "trace", "jac", "method", "hits", "misses", "mode"]
field(line, k) = match(Regex("\"$k\":(\"[^\"]*\"|[^,}]*)"), line).captures[1]
rows = [Dict(k => field(line, k) for k in keys_) for line in split(txt, "\n") if startswith(line, "{") || startswith(line, "[{")]
f(s) = tryparse(Float64, s)
ok = filter(r -> r["status"] == "\"ok\"", rows)
sp = [f(r["speed"]) for r in ok if f(r["speed"]) !== nothing]
println("cases = ", length(rows), "  ok = ", length(ok), "  with a speed = ", length(sp))
println("speedup: min ", minimum(sp), "  median ", median(sp), "  geomean ", exp(mean(log.(sp))), "  max ", maximum(sp))
println("speedup ranges: <1: ", count(<(1), sp), "  1-3: ", count(x -> 1 <= x < 3, sp), "  3-10: ", count(x -> 3 <= x < 10, sp),
    "  10-30: ", count(x -> 10 <= x < 30, sp), "  30-100: ", count(x -> 30 <= x < 100, sp), "  >=100: ", count(>=(100), sp))
er = [f(r["err"]) for r in ok if f(r["err"]) !== nothing]
println("error: max ", maximum(er), "  median ", median(er), "  above 1e-10: ", count(>(1e-10), er))
fs = [f(r["first"]) for r in ok]
println("first call, s: median ", median(fs), "  max ", maximum(fs), "  <1: ", count(<(1), fs), "  1-5: ", count(x -> 1 <= x < 5, fs), "  >=5: ", count(>=(5), fs))
h = [f(r["hit"]) for r in ok if f(r["hit"]) !== nothing]
println("cached call, us: median ", median(h), "  min ", minimum(h), "  max ", maximum(h))
println("fastest:"); for r in sort(ok, by = r -> -(f(r["speed"]) === nothing ? 0.0 : f(r["speed"])))[1:8]; println("  ", r["name"], "  ", r["speed"], "x"); end
println("slowest:"); for r in sort([r for r in ok if f(r["speed"]) !== nothing], by = r -> f(r["speed"]))[1:10]; println("  ", r["name"], "  ", r["speed"], "x  guards ", r["guards"], "  ", r["mode"], "  ", r["method"]); end
pk = Dict{String,Vector{String}}()
for r in rows; push!(get!(pk, r["pkg"], String[]), r["status"]); end
for (p, v) in sort(collect(pk)); println("  ", p, ": ", count(==("\"ok\""), v), "/", length(v)); end
println("reused a trace on a perturbed input: ", count(r -> f(r["hits"]) > 0, ok), "   recorded a new one on a perturbed input: ", count(r -> f(r["misses"]) > 0, ok))
println("guards: median ", median([f(r["guards"]) for r in ok]), "  max ", maximum([f(r["guards"]) for r in ok]), "   largest Jacobian nodes ", maximum([f(r["jac"]) for r in ok]))

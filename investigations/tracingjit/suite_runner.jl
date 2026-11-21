# Test suite runner for the tracing prototype: many functions from many registry packages.
#   julia --project=<env> suite_runner.jl list                 print the case names
#   julia --project=<env> suite_runner.jl run <first> <last>   run cases first..last (by position in the list)
#   julia --project=<env> suite_runner.jl run name1 name2 ...  run these cases
# One result line per case, also appended to the file named by the environment variable SUITE_RESULTS.
#
# A case "works" if the plain function runs, tracing it succeeds, and its Jacobian agrees with ForwardDiff (or, where
# ForwardDiff cannot differentiate the function, with finite differences).

include(joinpath(@__DIR__, "prototype.jl"))
using .TracingPrototype
const TP = TracingPrototype
using ForwardDiff, FiniteDifferences
using LinearAlgebra, Random, Printf, Statistics

# every package is imported, not used, so that the cases write qualified names and nothing clashes
const PACKAGES = ["SpecialFunctions", "NaNMath", "Distances", "StatsBase", "Polynomials", "DataInterpolations", "Interpolations",
    "BSplineKit", "QuadGK", "HCubature", "Roots", "IterativeSolvers", "Krylov", "Optim", "LsqFit", "OrdinaryDiffEqTsit5", "TaylorIntegration",
    "Distributions", "Loess", "Bessels", "Geodesy", "CoordinateTransformations", "Quaternions", "Rotations", "StaticArrays", "PDMats",
    "LogExpFunctions", "StatsFuns", "GenericLinearAlgebra", "HypergeometricFunctions", "SimpleNonlinearSolve", "FastGaussQuadrature"]
for p in PACKAGES
    try
        Core.eval(Main, :(import $(Symbol(p))))
    catch e
        println(stderr, "could not load ", p, ": ", first(split(sprint(showerror, e), '\n')))
    end
end

# Tracked methods for SpecialFunctions that FD has as Node operations
if isdefined(Main, :SpecialFunctions)
    for f in (:erf, :erfc, :gamma, :loggamma, :digamma)
        Core.eval(Main, :(SpecialFunctions.$f(a::$(TP.Tracked)) = $(TP.Tracked)(SpecialFunctions.$f(a.val), SpecialFunctions.$f(a.node))))
    end
end

struct Case
    name::String
    pkg::String
    f::Any
    args::Any            # a function returning the tuple of arguments
    wrt::Vector{Int}     # which arguments to differentiate (all if empty)
    extra::Vector{Any}   # more argument functions, to check other branches
    expect_failure::Bool # known not to work for a structural reason (for example LAPACK only)
end
const CASES = Case[]
add!(name, pkg, f, mk; wrt=Int[], extra=Any[], expect_failure=false) =
    push!(CASES, Case(name, pkg, f, () -> mk(rng_for(name)), wrt, Any[(() -> m(rng_for(name * "x$i"))) for (i, m) in enumerate(extra)], expect_failure))

# helpers for building arguments. Every case gets its own deterministic random generator.
rng_for(name) = MersenneTwister(hash(name) % UInt32)

for file in ("suite_cases_a.jl", "suite_cases_b.jl", "suite_cases_c.jl")
    path = joinpath(@__DIR__, file)
    isfile(path) && include(path)
end

# ---------------------------------------------------------------------------------------------------- evaluating a case
function flat_function(f, args, wrt)
    sizes = [a isa Number ? nothing : size(a) for a in args]
    ranges = TP.arg_ranges(args)
    return function (z)
        pos = 0
        newargs = ntuple(length(args)) do i
            if i in wrt
                len = args[i] isa Number ? 1 : length(args[i])
                piece = z[(pos+1):(pos+len)]
                pos += len
                args[i] isa Number ? piece[1] : reshape(piece, sizes[i])
            else
                args[i]
            end
        end
        out = f(newargs...)
        return out isa Number ? [out] : collect(vec(out))
    end
end
wrt_vector(args, wrt) = Float64[v for i in wrt for v in (args[i] isa Number ? (args[i],) : vec(args[i]))]
oneline(e) = (s = first(split(sprint(showerror, e), '\n')); length(s) > 200 ? s[1:200] * "..." : s)

function reference(f, args, wrt)
    g = flat_function(f, args, wrt)
    x0 = wrt_vector(args, wrt)
    try
        return ForwardDiff.jacobian(g, x0), :forwarddiff
    catch e1
        try
            return FiniteDifferences.jacobian(central_fdm(5, 1), g, x0)[1], :finitediff
        catch e2
            return nothing, :none
        end
    end
end

function timeit(thunk; budget=0.2)
    thunk(); thunk()
    n = 1
    while true
        t = @elapsed for _ in 1:n
            thunk()
        end
        t > budget && return t / n
        n *= 4
        n > 1_000_000 && return t / (n ÷ 4)
    end
end

Base.@kwdef mutable struct Result
    name::String = ""
    pkg::String = ""
    status::Symbol = :unknown   # :ok :mismatch :plain_failed :trace_failed :too_big :args_failed :no_output
    message::String = ""
    error::Float64 = NaN
    reference::Symbol = :none
    guards::Int = 0
    trace_nodes::Int = 0
    jacobian_nodes::Int = 0
    method::Symbol = :none
    first_call_s::Float64 = NaN
    t_hit_us::Float64 = NaN
    t_fd_us::Float64 = NaN
    t_plain_us::Float64 = NaN
    hits::Int = 0
    misses::Int = 0
    mode::Symbol = :none
    extra_ok::String = ""
end

function evaluate_case(c::Case)
    r = Result(name=c.name, pkg=c.pkg)
    args = try
        c.args()
    catch e
        r.status = :args_failed
        r.message = oneline(e)
        return r
    end
    wrt = isempty(c.wrt) ? collect(1:length(args)) : c.wrt
    out = try
        c.f(args...)
    catch e
        r.status = :plain_failed
        r.message = oneline(e)
        return r
    end
    if !(out isa Real || (out isa AbstractArray && all(x -> x isa Real, out)))
        r.status = :no_output
        r.message = "the function returned a $(typeof(out)), not a number or an array of numbers"
        return r
    end
    nout = out isa Number ? 1 : length(out)
    nwrt = sum(a -> a isa Number ? 1 : length(a), args[wrt])
    ref, how = reference(c.f, args, wrt)
    r.reference = how

    tj = TracingJacobian(c.f; wrt=wrt, in_place=true, max_traces=300)
    result = zeros(nout, nwrt)
    t = @elapsed try
        tj(result, args...)
    catch e
        r.status = occursin("TRACE_TOO_BIG", sprint(showerror, e)) ? :too_big : :trace_failed
        r.message = oneline(e)
        return r
    end
    r.first_call_s = t
    tr = tj.traces[1]
    r.guards = tr.stats.guards
    r.trace_nodes = tr.stats.trace_nodes
    r.jacobian_nodes = tr.stats.jacobian_nodes
    r.method = tr.stats.jacobian_method
    if ref !== nothing
        r.error = maximum(abs.(result .- ref) ./ max.(abs.(ref), 1.0))
        tol = how === :forwarddiff ? 1e-6 : 1e-3
        if !(r.error <= tol)
            r.status = :mismatch
            r.message = "max relative error $(round(r.error, sigdigits=2)) against $(how)"
        else
            r.status = :ok
        end
    else
        r.status = :ok
        r.message = "no reference derivative available"
    end

    # nearby inputs reuse the trace?
    h0, m0 = tj.hits, tj.misses
    rng = rng_for(c.name * "perturb")
    for _ in 1:8
        a2 = map(a -> a isa Number ? a * (1 + 1e-3 * randn(rng)) : a .* (1 .+ 1e-3 .* randn(rng, size(a)...)), args)
        try
            tj(result, a2...)
        catch
            break
        end
    end
    r.hits = tj.hits - h0
    r.misses = tj.misses - m0
    r.mode = tj.active

    # other branches
    notes = String[]
    for mk in c.extra
        try
            a2 = mk()
            ref2, how2 = reference(c.f, a2, wrt)
            tj(result, a2...)
            e2 = ref2 === nothing ? NaN : maximum(abs.(result .- ref2) ./ max.(abs.(ref2), 1.0))
            push!(notes, e2 <= (how2 === :forwarddiff ? 1e-6 : 1e-3) ? "ok" : "ERR$(round(e2, sigdigits=2))")
        catch e
            push!(notes, "failed")
        end
    end
    r.extra_ok = join(notes, ",")

    # speed
    try
        x = TP.flatten(args)
        r.t_hit_us = 1e6 * timeit(() -> TP.eval_jacobian!(tj, result, x))
        r.t_plain_us = 1e6 * timeit(() -> c.f(args...))
        if ref !== nothing && how === :forwarddiff
            g = flat_function(c.f, args, wrt)
            xw = wrt_vector(args, wrt)
            y = g(xw)
            g!(o, z) = (o .= g(z); nothing)
            cfg = ForwardDiff.JacobianConfig(g!, y, xw)
            Jf = similar(result)
            r.t_fd_us = 1e6 * timeit(() -> ForwardDiff.jacobian!(Jf, g!, y, xw, cfg))
        end
    catch e
        r.message *= " (timing failed: $(oneline(e)))"
    end
    return r
end

function format_result(r::Result)
    speed = isfinite(r.t_hit_us) && isfinite(r.t_fd_us) ? @sprintf("%.2fx", r.t_fd_us / r.t_hit_us) : "-"
    return @sprintf("%-8s %-34s %-18s err %-8s ref %-11s guards %5d trace %6d jac %7d [%-9s] first %7.2fs hit %9.2fus FD %9.2fus (%s) hits/misses %d/%d %s extra[%s] %s",
        string(r.status), r.name, r.pkg, isnan(r.error) ? "-" : @sprintf("%.1e", r.error), string(r.reference), r.guards, r.trace_nodes, r.jacobian_nodes,
        string(r.method), isnan(r.first_call_s) ? 0.0 : r.first_call_s, r.t_hit_us, r.t_fd_us, speed, r.hits, r.misses, string(r.mode), r.extra_ok, r.message)
end

function main()
    cmd = isempty(ARGS) ? "list" : ARGS[1]
    if cmd == "list"
        for (i, c) in enumerate(CASES)
            println(i, "\t", c.pkg, "\t", c.name)
        end
        return
    end
    selected = if length(ARGS) == 3 && all(a -> tryparse(Int, a) !== nothing, ARGS[2:3])
        CASES[parse(Int, ARGS[2]):min(parse(Int, ARGS[3]), length(CASES))]
    else
        filter(c -> c.name in ARGS[2:end], CASES)
    end
    TP.MAX_STATEMENTS[] = 400_000
    out = get(ENV, "SUITE_RESULTS", "")
    for c in selected
        line = try
            format_result(evaluate_case(c))
        catch e
            @sprintf("%-8s %-34s %-18s %s", "crashed", c.name, c.pkg, oneline(e))
        end
        println(line)
        flush(stdout)
        isempty(out) || open(out, "a") do io
            println(io, line)
        end
    end
end
main()

# Random expression graphs for checking `jacobian` against exact derivatives.
#
# A graph is described by a `GraphSpec`: `nvars` variables numbered 1:nvars,
# then one node per operation, numbered nvars+1, nvars+2, ...  Operations are
# +, -, *, sin, cos and squaring. Each operation's first operand is one of the
# few most recent nodes, so graphs are deep and share subexpressions heavily,
# and some roots are interior nodes of other roots. These are the conditions
# under which factoring has to handle overlapping subgraphs.
#
# The reference derivatives are computed by forward-mode propagation over the
# spec itself, so they do not depend on any FastDifferentiation code.
# `check_spec` builds the graph with FastDifferentiation, calls `jacobian`,
# evaluates the result at random points and compares it with the reference.
# Graphs whose exact values overflow Float64 are skipped, not counted as
# failures. `spec_code` prints Julia code that rebuilds a failing graph, for a
# bug report or a regression test.
#
# How it is used for testing
# --------------------------
# The "random graphs" test item in test/test-random-graphs.jl runs 3,220
# graphs from a fixed seed as part of the normal test suite (2000 tiny,
# 1000 small, 200 medium, 20 large; sizes as in the table below). The graphs
# take under a second; most of the item's run time is compilation.
#
# Running it on larger graphs
# ---------------------------
# If you suspect a problem, for example after changing factoring or code
# generation, run the fuzzer for longer, on larger graphs, with new seeds:
#
#     julia --project=. scripts/fuzz_factoring.jl [seed] [seconds_per_size]
#
# It runs each size in its SIZES table for `seconds_per_size` seconds
# (default 20), prints the time per graph and the number of failures, and
# prints reproduction code for the first few failures of each size. Edit
# SIZES there to try other shapes. To check one graph interactively:
#
#     include("test/random_graphs.jl"); using .RandomGraphs, Random
#     spec = random_spec(MersenneTwister(7); nvars=5, nops=50, nroots=4)
#     check_spec(spec)          # (:ok, "") or (:wrong / :error, message)
#     println(spec_code(spec))
#
# Timings and failure rates
# -------------------------
# Measured on 2026-11-14 with Julia 1.13 on a Windows laptop: 15 s per size
# for this branch and 10 s per size for `main` (before the factoring fix of
# PR #129). The branch is faster per graph because `main` falls into slow
# fallback paths on graphs it gets wrong.
#
#   size     vars   ops      roots   ms/graph       main failure rate      this branch
#                                    branch  main                          failures
#   tiny     2-3    3-6      1-3       0.09   0.13    4%  (3178/75338)     0 of 160,833
#   small    2-5    8-15     1-4       0.23   0.47   32%  (6743/21204)     0 of  64,371
#   medium   3-8    20-40    2-6       1.0   10.7    88%  (841/960)        0 of  15,316
#   large    4-12   60-120   3-10      9     305    100%  (33/33)          0 of   1,468
#   x-large  8-20   200-400  5-20    156   12093    100%  (2/2)            0 of      97
#
# A few thousand tiny and small graphs are enough to catch a bug of the kind
# fixed in PR #129, which is why the CI test uses mostly those sizes.

module RandomGraphs

import FastDifferentiation as FD
using Random

export GraphSpec, random_spec, check_spec, spec_code, fuzz

const BINARY = (:+, :-, :*)
const UNARY = (:sin, :cos, :sq)

struct GraphSpec
    nvars::Int
    ops::Vector{Tuple{Symbol,Int,Int}}   # (op, a, b); b == 0 for unary ops
    roots::Vector{Int}                   # node indices used as roots, in order
end

"""
    random_spec(rng; nvars, nops, nroots, unary_fraction=0.2, locality=6)

Random expression graph with `nops` operations on `nvars` variables. The first
operand of each operation is one of the `locality` most recent nodes, which
makes the graph deep and makes later nodes share earlier ones; the second is
uniform over all earlier nodes. Roots are the last node plus `nroots - 1`
other distinct operation nodes, so some roots are also interior nodes of other
roots."""
function random_spec(rng::AbstractRNG; nvars::Int, nops::Int, nroots::Int, unary_fraction=0.2, locality=6)
    ops = Tuple{Symbol,Int,Int}[]
    for i in 1:nops
        m = nvars + i - 1
        a = max(1, m - rand(rng, 0:locality-1))
        if rand(rng) < unary_fraction
            push!(ops, (rand(rng, UNARY), a, 0))
        else
            b = rand(rng, 1:m)
            b == a && (b = a == 1 ? min(2, m) : a - 1)
            push!(ops, (rand(rng, BINARY), a, b))
        end
    end
    last_node = nvars + nops
    others = shuffle(rng, collect(nvars+1:last_node-1))
    roots = vcat([last_node], others[1:min(nroots - 1, length(others))])
    return GraphSpec(nvars, ops, shuffle(rng, roots))
end

apply(op::Symbol, a, b) =
    op === :+ ? a + b :
    op === :- ? a - b :
    op === :* ? a * b :
    op === :sin ? sin(a) :
    op === :cos ? cos(a) :
    op === :sq ? a * a : error("unknown op $op")

"FastDifferentiation nodes for every node of `spec`."
function build(spec::GraphSpec, vars)
    nodes = FD.Node[vars...]
    for (op, a, b) in spec.ops
        push!(nodes, apply(op, nodes[a], b == 0 ? nothing : nodes[b]))
    end
    return nodes
end

"Exact Jacobian of the roots of `spec` at `x`, by forward-mode propagation."
function reference_jacobian(spec::GraphSpec, x::Vector{Float64})
    J = zeros(length(spec.roots), spec.nvars)
    for j in 1:spec.nvars
        v = copy(x)
        d = [i == j ? 1.0 : 0.0 for i in 1:spec.nvars]
        for (op, a, b) in spec.ops
            if op === :+
                push!(v, v[a] + v[b]); push!(d, d[a] + d[b])
            elseif op === :-
                push!(v, v[a] - v[b]); push!(d, d[a] - d[b])
            elseif op === :*
                push!(v, v[a] * v[b]); push!(d, d[a] * v[b] + v[a] * d[b])
            elseif op === :sin
                push!(v, sin(v[a])); push!(d, cos(v[a]) * d[a])
            elseif op === :cos
                push!(v, cos(v[a])); push!(d, -sin(v[a]) * d[a])
            else # :sq
                push!(v, v[a]^2); push!(d, 2 * v[a] * d[a])
            end
        end
        for (k, root) in enumerate(spec.roots)
            J[k, j] = d[root]
        end
    end
    return J
end

"Numeric value of a FastDifferentiation expression; `env` maps variable symbols to values."
function evaluate(n::FD.Node, env::Dict{Symbol,Float64}, memo::IdDict{FD.Node,Float64})
    cached = get(memo, n, nothing)
    cached === nothing || return cached
    result = if FD.is_variable(n)
        env[FD.value(n)]
    elseif FD.is_constant(n)
        Float64(FD.value(n))
    elseif FD.is_NoOp(n)
        evaluate(FD.children(n)[1], env, memo)
    else
        FD.value(n)(map(c -> evaluate(c, env, memo), FD.children(n))...)
    end
    memo[n] = result
    return result
end

"""
    check_spec(spec; rng, npoints=2, rtol=1e-8)

Returns `(:ok, "")`, `(:wrong, message)` if the Jacobian differs from the
reference at some point, `(:error, message)` if `jacobian` throws, or
`(:skip, "")` if the exact derivatives overflow `Float64` (long chains of
products can), in which case nothing can be compared."""
function check_spec(spec::GraphSpec; rng::AbstractRNG=Random.default_rng(), npoints=2, rtol=1e-8)
    points = [0.3 .+ 0.7 .* randn(rng, spec.nvars) for _ in 1:npoints]
    references = try
        [reference_jacobian(spec, x) for x in points]
    catch e
        e isa DomainError || rethrow()
        return (:skip, "")  # sin or cos of an overflowed value
    end
    all(R -> all(v -> isfinite(v) && abs(v) < 1e100, R), references) || return (:skip, "")

    vars = FD.make_variables(:v, spec.nvars)
    nodes = build(spec, vars)
    jac = try
        FD.jacobian(nodes[spec.roots], vars)
    catch e
        return (:error, sprint(showerror, e))
    end
    for (x, R) in zip(points, references)
        env = Dict(FD.value(vars[i]) => x[i] for i in 1:spec.nvars)
        memo = IdDict{FD.Node,Float64}()
        J = try
            map(n -> evaluate(n, env, memo), jac)
        catch e
            return (:wrong, "evaluating the Jacobian at x = $x threw $(sprint(showerror, e))")
        end
        scale = max(1.0, maximum(abs, R; init=0.0))
        if maximum(abs.(J .- R); init=0.0) > rtol * scale
            return (:wrong, "at x = $x\n got      $J\n expected $R")
        end
    end
    return (:ok, "")
end

"Julia source that rebuilds `spec`, for pasting into a bug report or a regression test."
function spec_code(spec::GraphSpec)
    name(i) = i <= spec.nvars ? "v$i" : "n$i"
    io = IOBuffer()
    println(io, "@variables ", join(["v$i" for i in 1:spec.nvars], " "))
    for (k, (op, a, b)) in enumerate(spec.ops)
        lhs = name(spec.nvars + k)
        rhs = op === :sq ? "$(name(a)) * $(name(a))" :
              b == 0 ? "$(op)($(name(a)))" : "$(name(a)) $(op) $(name(b))"
        println(io, lhs, " = ", rhs)
    end
    print(io, "jacobian([", join(name.(spec.roots), ", "), "], [", join(["v$i" for i in 1:spec.nvars], ", "), "])")
    return String(take!(io))
end

"""
    fuzz(rng, ntrials; nvars, nops, nroots, kwargs...)

Checks `ntrials` random graphs. Returns the failures as `(spec, status, message)`
tuples and the number of graphs skipped because their derivatives overflow."""
function fuzz(rng::AbstractRNG, ntrials::Integer; nvars, nops, nroots, kwargs...)
    failures = Tuple{GraphSpec,Symbol,String}[]
    skipped = 0
    for _ in 1:ntrials
        spec = random_spec(rng; nvars=rand(rng, nvars), nops=rand(rng, nops), nroots=rand(rng, nroots), kwargs...)
        status, msg = check_spec(spec; rng)
        if status === :skip
            skipped += 1
        elseif status !== :ok
            push!(failures, (spec, status, msg))
        end
    end
    return failures, skipped
end

end # module

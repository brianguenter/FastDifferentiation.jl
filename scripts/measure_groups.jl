#!/usr/bin/env julia
# Reproduces the measurements table in HOW_FASTDIFFERENTIATION_WORKS.md, section 14 (Cost).
#
# Usage (from the repository root):
#     julia --project=. scripts/measure_groups.jl
#
# For each problem it replays `factor!`, and just before each call to
# `factor_subgraph!` runs the group dynamic program (`subgraph_path_groups`) on
# the subgraph about to be factored. It reports the number of steps, the
# largest number of groups at the dominated node c, the largest table |G(u)|
# at any node, and the (n, m) = (|dom(S)|, |nondom(S)|) of the subgraph with
# the most groups.
#
# The random graphs use MersenneTwister(1). Julia does not guarantee that
# `rand` streams are identical across versions; the table was produced with
# Julia 1.13.

import FastDifferentiation as FD
using Random

function measure(name, roots)
    gr = FD.DerivativeGraph(roots)
    heap = FD.compute_factorable_subgraphs(gr)
    nsteps = 0; max_groups = 0; max_table = 0; worst = (0, 0, 0)
    while !isempty(heap)
        sub = pop!(heap)
        nsteps += 1
        region = FD.subgraph_region(sub)
        if sub isa FD.FactorableSubgraph{Int64,FD.DominatorSubgraph}
            memo = Dict{Int64,FD.PathGroups}()
            groups = FD.subgraph_path_groups(sub, FD.dominated_node(sub), memo)
            tbl = maximum(length.(values(memo)); init=0)
        else
            groups = FD.subgraph_path_groups(sub, region)
            tbl = length(groups)
        end
        g = length(groups)
        n = count(FD.reachable_dominance(sub)); m = count(FD.non_dominance_mask(sub))
        if g > max_groups
            max_groups = g; worst = (n, m, g)
        end
        max_table = max(max_table, tbl)
        FD.factor_subgraph!(sub)
    end
    println(rpad(name, 34), " steps=", lpad(nsteps, 6), "  max groups=", lpad(max_groups, 3),
            "  max |G(u)|=", lpad(max_table, 3), "  (n, m) at max groups = ", worst[1:2])
end

FD.@variables x y
p = x + y; q = x + p; r = p + q
measure("worked example", [q, r])

# Rosenbrock Hessian (graph of the gradient)
x = FD.make_variables(:x, 60)
ros = sum(100 * (x[i+1] - x[i]^2)^2 + (1 - x[i])^2 for i in 1:59)
measure("Rosenbrock gradient, n=60", vec(FD.jacobian([ros], x)))

# chain sum sin(x_i x_{i+1}) exp(x_i), Hessian
y = FD.make_variables(:y, 40)
ch = sum(sin(y[i] * y[i+1]) * exp(y[i]) for i in 1:39)
measure("chain gradient, n=40", vec(FD.jacobian([ch], y)))

# dense map x_i * sin(sum_j x_j (i+j)/n)
z = FD.make_variables(:z, 20)
measure("dense map, n=20", [z[i] * sin(sum(z[j] * (i + j) / 20 for j in 1:20)) for i in 1:20])

# random expression graphs with heavy sharing
function random_graph(rng, nin, nout, nops)
    vs = FD.make_variables(:r, nin)
    pool = Any[vs...]
    for _ in 1:nops
        a = pool[max(1, length(pool) - rand(rng, 0:5))]
        b = pool[rand(rng, 1:length(pool))]
        op = rand(rng, 1:6)
        push!(pool, op == 1 ? a + b : op == 2 ? a - b : op == 3 ? a * b : op == 4 ? sin(a) : op == 5 ? a / (1 + b^2) : a^2)
    end
    return [pool[end - k] for k in 0:nout-1]
end
rng = MersenneTwister(1)
for t in 1:4
    measure("random graph $t (8 in, 8 out, 70 ops)", random_graph(rng, 8, 8, 70))
end

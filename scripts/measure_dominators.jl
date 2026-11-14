#!/usr/bin/env julia
# Reproduces the dominator-table charts in HOW_FASTDIFFERENTIATION_WORKS.md,
# section 14 (Cost).
#
# Usage (from the repository root):
#     julia --project=. scripts/measure_dominators.jl          # writes img/fig-dom-*.svg
#     julia --project=. scripts/measure_dominators.jl outdir   # writes the charts elsewhere
#
# For each scalable example it builds the derivative graph at several sizes,
# runs the dominator and postdominator tables exactly as
# compute_factorable_subgraphs does (one per root, one per variable), and counts
# the steps taken by the `intersect` walks. The counter is added by redefining
# FastDifferentiation.intersection with the same loop plus a counter.
#
# It writes one chart per example (intersect steps against nodes + edges) and
# prints the measurements as a Markdown table.
#
# Examples: spherical harmonics, Chebyshev polynomials and the heat equation of
# issue #60 (all from test/runtests.jl); the Rosenbrock function, a chain of
# sin(x_i x_{i+1}) exp(x_i) terms and a dense map; and the linear
# model-predictive-control rollout of issue #26 (see scripts/scalable_examples.jl).

import FastDifferentiation as FD
include(joinpath(@__DIR__, "scalable_examples.jl"))
include(joinpath(@__DIR__, "growth_chart.jl"))

const STEPS = Ref(0)
@eval FD function intersection(order_test, node1::T, node2::T, idoms::Dict{T,T}, start_index::T) where {T<:Integer}
    while node1 != node2
        $STEPS[] += 1
        if order_test(node1, node2)
            node1 = idoms[node1]
        else
            node2 = idoms[node2]
        end
    end
    return node1
end


const EXAMPLES = [
    ("spherical harmonics", "fig-dom-spherical-harmonics.svg", s -> "order $s", spherical_harmonics, 4:4:24),
    ("Chebyshev polynomial", "fig-dom-chebyshev.svg", s -> "degree $s", chebyshev, 50:50:300),
    ("heat equation (issue #60)", "fig-dom-heat.svg", s -> "$s points", heat, 40:40:240),
    ("Rosenbrock function", "fig-dom-rosenbrock.svg", s -> "$s variables", rosenbrock, 20:20:120),
    ("Σ sin(xᵢxᵢ₊₁)·exp(xᵢ)", "fig-dom-sin-exp-chain.svg", s -> "$s variables", sin_exp_chain, 20:20:120),
    ("dense map xᵢ·sin(Σⱼ xⱼ(i+j)/n)", "fig-dom-dense-map.svg", s -> "n = $s", dense_map, 5:5:30),
    ("linear MPC rollout (issue #26)", "fig-dom-mpc.svg", s -> "horizon $s", mpc, 5:7:40),
]

# ---------------------------------------------------------------- measurement
function measure(roots)
    gr = FD.DerivativeGraph(roots)
    E = sum(length(FD.children(e)) for e in values(FD.edges(gr)))
    V = length(FD.nodes(gr))
    tables = FD.codomain_dimension(gr) + FD.domain_dimension(gr)
    STEPS[] = 0
    doms = Dict{Int,Int}()
    for i in 1:FD.codomain_dimension(gr)
        FD.compute_dom_table(gr, true, i, FD.root_index_to_postorder_number(gr, i), doms)
    end
    for j in 1:FD.domain_dimension(gr)
        FD.compute_dom_table(gr, false, j, FD.variable_index_to_postorder_number(gr, j), doms)
    end
    return (V=V, E=E, tables=tables, steps=STEPS[])
end


outdir = isempty(ARGS) ? joinpath(@__DIR__, "..", "img") : ARGS[1]
mkpath(outdir)
measure(chebyshev(5))   # compile before measuring

println("| example | size | nodes + edges | tables (roots + variables) | intersect steps | steps per (node + edge) per table |")
println("|:--|--:|--:|--:|--:|--:|")
for (title, file, sizephrase, f, sizes) in EXAMPLES
    data = [measure(f(s)) for s in sizes]
    for (s, d) in zip(sizes, data)
        println("| $title | $(sizephrase(s)) | $(d.V + d.E) | $(d.tables) | $(d.steps) | $(round(d.steps / ((d.V + d.E) * d.tables), digits=3)) |")
    end
    write(joinpath(outdir, file), growth_chart(title, [d.V + d.E for d in data], [d.steps for d in data];
        quantity="intersect steps", end_label="$(fmt(last(data).steps)) steps at $(sizephrase(last(sizes)))",
        zero_message="no intersect steps at any size"))
end

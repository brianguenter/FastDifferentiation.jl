#!/usr/bin/env julia
# Reproduces the group growth charts and the segment table in
# HOW_FASTDIFFERENTIATION_WORKS.md, Appendix A.
#
# Usage (from the repository root):
#     julia --project=. scripts/measure_segments.jl          # writes img/fig-groups-*.svg
#     julia --project=. scripts/measure_segments.jl outdir   # writes the charts elsewhere
#
# For each family of graphs and each size, it replays `factor!` and, just before
# each call to `factor_subgraph!`, runs the group tables (`subgraph_path_groups`)
# on the subgraph about to be factored. It records, over all steps of one graph:
#   - segments: the most segments any single step sees (the sum of n over its
#     groups, i.e. the chains a naive implementation would have to list),
#   - groups: the most groups at the dominated node of any step, and the
#     (n, m) of that step,
#   - table: the largest table |G(u)| at any node of any step.
# For random graphs it reports the worst graph of a fixed-seed sample at each
# size. It writes one chart per family (most groups against nodes + edges) and
# prints all measurements as a Markdown table.
#
# The largest random family (80 roots and 80 variables) takes several minutes.

include(joinpath(@__DIR__, "scalable_examples.jl"))
include(joinpath(@__DIR__, "growth_chart.jl"))
include(joinpath(@__DIR__, "..", "test", "random_graphs.jl"))
using .RandomGraphs, Random

function replay(roots)
    gr = FD.DerivativeGraph(roots)
    size = length(FD.nodes(gr)) + sum(length(FD.children(e)) for e in values(FD.edges(gr)))
    heap = FD.compute_factorable_subgraphs(gr)
    segs = 0; groups = 0; table = 0; nm = (0, 0)
    while !isempty(heap)
        sub = pop!(heap)
        region = FD.subgraph_region(sub)
        if sub isa FD.FactorableSubgraph{Int64,FD.DominatorSubgraph}
            memo = Dict{Int64,FD.PathGroups}()
            g = FD.subgraph_path_groups(sub, FD.dominated_node(sub), memo)
            table = max(table, maximum(length.(values(memo)); init=0))
        else
            g = FD.subgraph_path_groups(sub, region)
            table = max(table, length(g))
        end
        segs = max(segs, sum((x[4] for x in values(g)); init=0))
        if length(g) > groups
            groups = length(g)
            nm = (count(FD.reachable_dominance(sub)), count(FD.non_dominance_mask(sub)))
        end
        FD.factor_subgraph!(sub)
    end
    return (size=size, segs=segs, groups=groups, table=table, nm=nm)
end

"Worst of `trials` random graphs; `size` is the mean size of the graphs."
function random_worst(rng, trials; nvars, nops, nroots)
    rs = map(1:trials) do _
        spec = random_spec(rng; nvars=nvars(), nops=nops, nroots=nroots())
        replay(RandomGraphs.build(spec, FD.make_variables(:v, spec.nvars))[spec.roots])
    end
    worst = rs[argmax([r.groups for r in rs])]
    return (size=round(Int, sum(r.size for r in rs) / trials), segs=maximum(r.segs for r in rs),
            groups=worst.groups, table=maximum(r.table for r in rs), nm=worst.nm)
end

const FAMILIES = [
    ("k × k grid", "fig-groups-grid.svg", k -> "k = $k", k -> replay(grid(k)), 3:2:13),
    ("spherical harmonics", "fig-groups-spherical-harmonics.svg", s -> "order $s", s -> replay(spherical_harmonics(s)), 4:4:24),
    ("Chebyshev polynomial", "fig-groups-chebyshev.svg", s -> "degree $s", s -> replay(chebyshev(s)), 50:50:300),
    ("heat equation (issue #60)", "fig-groups-heat.svg", s -> "$s points", s -> replay(heat(s)), 40:40:240),
    ("Rosenbrock function", "fig-groups-rosenbrock.svg", s -> "$s variables", s -> replay(rosenbrock(s)), 20:20:120),
    ("Σ sin(xᵢxᵢ₊₁)·exp(xᵢ)", "fig-groups-sin-exp-chain.svg", s -> "$s variables", s -> replay(sin_exp_chain(s)), 20:20:120),
    ("dense map xᵢ·sin(Σⱼ xⱼ(i+j)/n)", "fig-groups-dense-map.svg", s -> "n = $s", s -> replay(dense_map(s)), 5:5:30),
    ("linear MPC rollout (issue #26)", "fig-groups-mpc.svg", s -> "horizon $s", s -> replay(mpc(s)), 5:7:40),
    ("random graphs, 3–8 roots and variables", "fig-groups-random.svg", s -> "$s operations",
        let rng = MersenneTwister(11)
            s -> random_worst(rng, 15; nvars=() -> rand(rng, 3:8), nops=s, nroots=() -> rand(rng, 3:8))
        end, [40, 160, 640, 2560]),
    ("random graphs, k roots and k variables", "fig-groups-random-wide.svg", k -> "k = $k",
        let rng = MersenneTwister(21)
            k -> random_worst(rng, k <= 40 ? 4 : 2; nvars=() -> k, nops=20k, nroots=() -> k)
        end, [5, 10, 20, 40, 80]),
]

outdir = isempty(ARGS) ? joinpath(@__DIR__, "..", "img") : ARGS[1]
mkpath(outdir)
replay(grid(3))   # compile before measuring

println("| graphs | size | nodes + edges | most segments in one step | most groups in one step | (n, m) of that step | largest table |")
println("|:--|:--|--:|--:|--:|:--|--:|")
for (title, file, sizephrase, measure, sizes) in FAMILIES
    rs = [measure(s) for s in sizes]
    for (s, r) in zip(sizes, rs)
        println("| $title | $(sizephrase(s)) | $(r.size) | $(r.segs) | $(r.groups) | $(r.nm) | $(r.table) |")
    end
    write(joinpath(outdir, file), growth_chart(title, [r.size for r in rs], [r.groups for r in rs];
        quantity="most groups in one step", end_label="$(last(rs).groups) group$(last(rs).groups == 1 ? "" : "s") at $(sizephrase(last(sizes)))", integer=true,
        zero_message="no factoring step had a group"))
end

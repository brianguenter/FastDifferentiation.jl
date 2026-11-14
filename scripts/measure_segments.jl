#!/usr/bin/env julia
# Reproduces the segment and group growth table in
# HOW_FASTDIFFERENTIATION_WORKS.md, Appendix A.
#
# Usage (from the repository root):
#     julia --project=. scripts/measure_segments.jl
#
# For each family of graphs and each size, it replays `factor!` and, just before
# each call to `factor_subgraph!`, runs the group tables (`subgraph_path_groups`)
# on the subgraph about to be factored. It records, over all steps of one graph:
#   - segments: the most segments any single step sees (the sum of n over its
#     groups, i.e. the chains a naive implementation would have to list),
#   - groups: the most groups at the dominated node of any step,
#   - table: the largest table |G(u)| at any node of any step.
# For random graphs it reports the worst graph out of a fixed-seed sample at
# each size. Each row lists one value per size, smallest size first.

include(joinpath(@__DIR__, "scalable_examples.jl"))
include(joinpath(@__DIR__, "..", "test", "random_graphs.jl"))
using .RandomGraphs, Random

function replay(roots)
    gr = FD.DerivativeGraph(roots)
    size = length(FD.nodes(gr)) + sum(length(FD.children(e)) for e in values(FD.edges(gr)))
    heap = FD.compute_factorable_subgraphs(gr)
    segs = 0; groups = 0; table = 0
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
        groups = max(groups, length(g))
        FD.factor_subgraph!(sub)
    end
    return (size=size, segs=segs, groups=groups, table=table)
end

function random_worst(rng, nops, trials)
    worst = (size=0, segs=0, groups=0, table=0)
    for _ in 1:trials
        spec = random_spec(rng; nvars=rand(rng, 3:8), nops=nops, nroots=rand(rng, 3:8))
        vars = FD.make_variables(:v, spec.nvars)
        r = replay(RandomGraphs.build(spec, vars)[spec.roots])
        worst = (size=max(worst.size, r.size), segs=max(worst.segs, r.segs),
                 groups=max(worst.groups, r.groups), table=max(worst.table, r.table))
    end
    return worst
end

arrow(v) = join(v, " → ")

const FAMILIES = [
    ("k × k grid (paths: 6 → 70 → 924 → 12,870 → 184,756 → 2,704,156)", "k", grid, 3:2:13),
    ("spherical harmonics", "order", spherical_harmonics, 4:4:24),
    ("Chebyshev polynomial", "degree", chebyshev, 50:50:300),
    ("heat equation (issue #60)", "points", heat, 40:40:240),
    ("Rosenbrock function", "variables", rosenbrock, 20:20:120),
    ("Σ sin(xᵢxᵢ₊₁)·exp(xᵢ)", "variables", sin_exp_chain, 20:20:120),
    ("dense map xᵢ·sin(Σⱼ xⱼ(i+j)/n)", "n", dense_map, 5:5:30),
    ("linear MPC rollout (issue #26)", "horizon", mpc, 5:7:40),
]

replay(grid(3))   # compile before measuring
println("| graphs | sizes | nodes + edges | most segments in one step | most groups | largest table |")
println("|:--|:--|:--|:--|:--|:--|")
for (label, word, f, sizes) in FAMILIES
    rs = [replay(f(s)) for s in sizes]
    println("| $label | $word $(first(sizes))–$(last(sizes)) | $(first(rs).size)–$(last(rs).size) | ",
        arrow([r.segs for r in rs]), " | ", arrow([r.groups for r in rs]), " | ", arrow([r.table for r in rs]), " |")
end
rng = MersenneTwister(11)
nops = (40, 160, 640, 2560)
rs = [random_worst(rng, n, 15) for n in nops]
println("| random graphs, worst of 15 per size | operations $(first(nops))–$(last(nops)) | $(first(rs).size)–$(last(rs).size) | ",
    arrow([r.segs for r in rs]), " | ", arrow([r.groups for r in rs]), " | ", arrow([r.table for r in rs]), " |")

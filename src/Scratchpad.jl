#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    @variables x

    nx = Node(x)
    gr = DerivativeGraph((cos(nx) * cos(nx)) + nx)
    Vis.draw_dot(gr)
    Vis.draw_dot(gr)
    sub = FactorableSubgraph{Int64,DominatorSubgraph}(gr, 4, 1, BitVector([1]), BitVector([1]), BitVector([1]))

    edges_4_1 = collect(subgraph_edges(sub))

    sub = FactorableSubgraph{Int64,PostDominatorSubgraph}(gr, 1, 4, BitVector([1]), BitVector([1]), BitVector([1]))
    edges_1_4 = collect(subgraph_edges(sub))

    @assert count(x -> vertices(x) == (4, 3), edges_4_1) == 1
    @assert count(x -> vertices(x) == (4, 1), edges_4_1) == 1
    @assert count(x -> vertices(x) == (3, 2), edges_4_1) == 2
    @assert count(x -> vertices(x) == (2, 1), edges_4_1) == 1

    @assert count(x -> vertices(x) == (4, 3), edges_1_4) == 1
    @assert count(x -> vertices(x) == (4, 1), edges_1_4) == 1
    @assert count(x -> vertices(x) == (3, 2), edges_1_4) == 2
    @assert count(x -> vertices(x) == (2, 1), edges_1_4) == 1
end
export test

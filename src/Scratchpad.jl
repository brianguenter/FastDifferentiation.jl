#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    @variables x

    nx = Node(x)
    gr = DerivativeGraph((cos(nx) * cos(nx)) + x)
    Vis.draw_dot(gr)
    sub = FactorableSubgraph{Int64,DominatorSubgraph}(gr, 4, 1, BitVector([1]), BitVector([1]), BitVector([1]))

    sub_edges = subgraph_edges(sub)
    computed = vertices.(sub_edges)
    correct = ((4, 1), ())
end
export test

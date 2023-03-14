#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using Profile
using PProf
using .TestCases

function test()
    doms = dominator_subgraph.((
        (graph, 4, 2, BitVector([1, 0]), BitVector([1, 0])BitVector([1])),
        (graph, 4, 1, BitVector([1, 1]), BitVector([1, 0])BitVector([1]))))
    pdoms = postdominator_subgraph.((
        (graph, 1, 3, BitVector([1, 1]), BitVector([1]), BitVector([1, 1])),
        (graph, 1, 4, BitVector([1, 1]), BitVector([1]), BitVector([1, 0]))))
    subs2 = collect((pdoms..., doms...))
end
export test

function profile()
    graph, qx, qy, qz = to_graph(25)
    Profile.clear()
    @profile symbolic_jacobian!(graph, [Node(qx), Node(qy), Node(qz)])
    pprof()
end
export profile

#changed
#change
export test

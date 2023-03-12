#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using .TestCases
using StaticArrays
using .SphericalHarmonics
using FiniteDifferences
using Profile
using PProf

function test()
    graph, four_2_subgraph, one_3_subgraph, _ = simple_dominator_dgraph()
end
export test

function profile()
    graph, qx, qy, qz = to_graph(25)
    Profile.clear()
    @profile symbolic_jacobian!(graph, [Node(qx), Node(qy), Node(qz)])
    pprof()
end
export profile

export test

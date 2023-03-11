#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using .TestCases
using StaticArrays
using .SphericalHarmonics
using FiniteDifferences
using Profile
using PProf

function test()
    @variables x y

    nv1 = Node(x)
    nv2 = Node(y)
    n3 = nv1 * nv2
    n4 = n3 * nv1
    n5 = n3 * n4

    graph = RnToRmGraph([n4, n5])
    # factor_subgraph!(graph, postdominator_subgraph(2, 4, 2, BitVector([0, 1]), BitVector([0, 1])))
    subs, _ = compute_factorable_subgraphs(graph)
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

#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()

    @variables x y

    nx1 = Node(x)
    ny2 = Node(y)
    n3 = nx1 * ny2
    n4 = n3 * ny2
    n5 = n3 * n4

    graph = DerivativeGraph([n5, n4])
    Vis.draw_dot(graph)
    readline()
    tmp = postdominator_subgraph(graph, 2, 4, BitVector([0, 1]), BitVector([0, 1]), BitVector([0, 1]))
    factor_subgraph!(tmp)
    Vis.draw_dot(graph)
end
export test


function make_subdata()
    @variables x y

    nx = Node(x)
    ny = Node(y)
    n2 = nx * ny
    n4 = n2 * ny
    n5 = n2 * n4


    graph = DerivativeGraph([n4, n5])
    subs, _ = compute_factorable_subgraphs(graph)

    _5_3 = subs[1]
    @assert (5, 3) == vertices(_5_3)
    e5_3 = edges(graph, 5, 3)[1]
    e3_4 = edges(graph, 3, 4)[1]

    return path(_5_3, e3_4)
end
export make_subdata

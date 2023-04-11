#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()

    @variables x y

    nv1 = Node(x)
    nv2 = Node(y)
    n3 = nv1 * nv2
    n4 = n3 * nv1
    n5 = n3 * n4

    graph = DerivativeGraph([n4, n5])
    # factor_subgraph!(graph, postdominator_subgraph(2, 4, 2, BitVector([0, 1]), BitVector([0, 1])))
    subs, subgraph_dict = compute_factorable_subgraphs(graph)

    _5_3 = dominator_subgraph(graph, 5, 3, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _1_4 = postdominator_subgraph(graph, 1, 4, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _3_5 = postdominator_subgraph(graph, 3, 5, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _4_1 = dominator_subgraph(graph, 4, 1, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _5_1 = dominator_subgraph(graph, 5, 1, Bool[0, 1], Bool[0, 1], Bool[1, 0])
    _1_5 = postdominator_subgraph(graph, 1, 5, Bool[1, 0], Bool[0, 1], Bool[1, 0])

    sub_eval = evaluate_subgraph(_5_3)
    factor_subgraph!(_5_3)
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

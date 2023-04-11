#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
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

    add_non_dom_edges!(_5_3)
    #single edge 3,4 should be split into two: ([r1,r2],[v1,v2]) -> ([r1],[v1,v2]),([r2],[v1,v2])
    edges3_4 = edges(graph, 4, 3)
    @assert length(edges3_4) == 2
    test_edge = PathEdge(4, 3, ny, BitVector([1, 1]), BitVector([0, 1]))
    @assert count(value_equal.(edges3_4, Ref(test_edge))) == 1
    test_edge = (PathEdge(4, 3, ny, BitVector([1, 1]), BitVector([1, 0])))
    @assert count(value_equal.(edges3_4, Ref(test_edge))) == 1

    graph = DerivativeGraph([n4, n5])
    subs, _ = compute_factorable_subgraphs(graph)
    _2_4 = subs[2]
    @assert (2, 4) == vertices(_2_4)

    add_non_dom_edges!(_2_4)
    #single edge 3,4 should be split in two: ([r1,r2],[v1,v2])->([r1,r2],[v1]),([r1,r2],[v2])
    edges3_4 = edges(graph, 4, 3)
    @assert length(edges3_4) == 2
    test_edge = PathEdge(4, 3, ny, BitVector([1, 0]), BitVector([1, 1]))
    @assert count(value_equal.(edges3_4, Ref(test_edge))) == 1
    test_edge = (PathEdge(4, 3, ny, BitVector([0, 1]), BitVector([1, 1])))
    @assert count(value_equal.(edges3_4, Ref(test_edge))) == 1
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

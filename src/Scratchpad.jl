#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    function edge_fields_equal(edge1, edge2)
        return edge1.top_vertex == edge2.top_vertex &&
               edge1.bott_vertex == edge2.bott_vertex &&
               edge1.edge_value == edge2.edge_value &&
               edge1.reachable_variables == edge2.reachable_variables &&
               edge1.reachable_roots == edge2.reachable_roots
    end

    @variables x y

    nx = Node(x)
    ny = Node(y)
    n2 = nx * ny
    n4 = n2 * ny
    n5 = n2 * n4

    graph = DerivativeGraph([n4, n5])
    subs_heap = compute_factorable_subgraphs(graph)
    subs = extract_all!(subs_heap)
    _5_3 = subs[1]
    @assert (5, 3) == vertices(_5_3)

    add_non_dom_edges!(_5_3)
    #single edge 3,4 should be split into two: ([r1,r2],[v1,v2]) -> ([r1],[v1,v2]),([r2],[v1,v2])
    edges3_4 = edges(graph, 4, 3)
    @assert length(edges3_4) == 2
    test_edge = PathEdge(4, 3, ny, BitVector([1, 1]), BitVector([0, 1]))
    @assert count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1
    test_edge = (PathEdge(4, 3, ny, BitVector([1, 1]), BitVector([1, 0])))
    @assert count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1

    graph = DerivativeGraph([n4, n5])
    sub_heap = compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)
    _2_4 = subs[2]
    @assert (2, 4) == vertices(_2_4)

    add_non_dom_edges!(_2_4)
    #single edge 3,4 should be split in two: ([r1,r2],[v1,v2])->([r1,r2],[v1]),([r1,r2],[v2])
    edges3_4 = edges(graph, 4, 3)
    @assert length(edges3_4) == 2
    test_edge = PathEdge(4, 3, ny, BitVector([1, 0]), BitVector([1, 1]))
    @assert count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1
    test_edge = (PathEdge(4, 3, ny, BitVector([0, 1]), BitVector([1, 1])))
    @assert count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1
    return nothing
end
export test

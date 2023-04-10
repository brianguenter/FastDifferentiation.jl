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

    subs = compute_factorable_subgraphs(graph)[1]
    _5_3 = subs[findfirst(x -> x.subgraph == (5, 3), subs)]
    _4_1 = subs[findfirst(x -> x.subgraph == (4, 1), subs)]
    _1_4 = subs[findfirst(x -> x.subgraph == (1, 4), subs)]

    _5_3_edges = ((5, 3), (4, 3), (5, 4))
    _4_1_edges = ((4, 1), (4, 3), (3, 1))
    _1_4_edges = ((4, 1), (4, 3), (3, 1))

    @assert length(edges(_5_3)) == 3
    @assert length(edges(_4_1)) == 3
    @assert length(edges(_1_4)) == 3

    function test_edges(subgraph, edge_list)
        for edge in edges(subgraph)
            @assert (top_vertex(edge), bott_vertex(edge)) in edge_list
        end
    end

    for (subgraph, edges_of_subgraph) in zip((_5_3, _4_1, _1_4), (_5_3_edges, _4_1_edges), (_1_4_edges))
        test_edges(subgraph, edges_of_subgraph)
    end
end
export test

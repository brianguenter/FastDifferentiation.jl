#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables x

    nx = Node(x)
    func = nx * nx

    graph = DerivativeGraph([func])
    subs, _ = compute_factorable_subgraphs(graph)

    test_sub = subs[1]
    edges = parent_edges(graph, dominated_node(test_sub))
    rroots = reachable_roots(edges[1])
    rroots .= rroots .& .!rroots

    @assert !connected_path(test_sub, edges[1])
    @assert connected_path(test_sub, edges[2])
end
export test

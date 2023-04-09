#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables x y

    nx = Node(x)
    func = nx * nx

    gr = DerivativeGraph([func])
    subs, _ = compute_factorable_subgraphs(gr)

    test_sub = subs[1]
    etmp = parent_edges(gr, dominated_node(test_sub))
    rroots = reachable_roots(etmp[1])
    rroots .= rroots .& .!rroots

    @assert !connected_path(test_sub, etmp[1])
    @assert connected_path(test_sub, etmp[2])
end
export test

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
    subs_heap = compute_factorable_subgraphs(gr)
end
export test

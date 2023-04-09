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
    _2_4 = subs[2]
    _3_5 = subs[4]

    etmp = edges(graph, 3, 5)[1]
    @assert connected_path(_5_3, etmp)


    etmp = edges(graph, 3, 4)[1]
    @assert connected_path(_5_3, etmp)
    rts = reachable_roots(etmp)
    rts[2] = 0

    @assert !connected_path(_5_3, etmp)
    #reset path
    rts[2] = 1

    e2_4 = edges(graph, 2, 4)[1]
    @assert connected_path(_2_4, e2_4)
    e2_3 = edges(graph, 2, 3)[1]
    @assert connected_path(_2_4, e2_3)
    e3_4 = edges(graph, 3, 4)[1]
    vars = reachable_variables(e3_4)
    @. vars &= !vars
    @assert !connected_path(_2_4, e3_4)
end
export test

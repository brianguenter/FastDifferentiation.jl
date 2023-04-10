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
    println(subs)
    _4_2 = subs[3]
    @assert (4, 2) == vertices(_4_2)

    add_non_dom_edges!(_4_2)
    Vis.draw_dot(graph)
    readline()

    graph = DerivativeGraph([n4, n5])
    subs, _ = compute_factorable_subgraphs(graph)
    _3_5 = subs[4]
    @assert (3, 5) == vertices(_3_5)

    add_non_dom_edges!(_3_5)
    Vis.draw_dot(graph)
    # _5_3 = subs[1]
    # @assert (5, 3) == vertices(_5_3) #these are not tests. Put these here in case some future change to code causes order of subgraphs to change. Shouldn't happen but could.
    # _2_4 = subs[2]
    # @assert (2, 4) == vertices(_2_4)
    # _3_5 = subs[4]
    # @assert (3, 5) == vertices(_3_5)



    # etmp = edges(graph, 3, 4)[1]
    # @test connected_path(_5_3, etmp)
    # rts = reachable_roots(etmp)
    # rts[2] = 0

    # @test !connected_path(_5_3, etmp)
    # #reset path
    # rts[2] = 1

    # e2_4 = edges(graph, 2, 4)[1]
    # @test connected_path(_2_4, e2_4)
    # e2_3 = edges(graph, 2, 3)[1]
    # @test connected_path(_2_4, e2_3)
    # e3_4 = edges(graph, 3, 4)[1]
    # vars = reachable_variables(e3_4)
    # @. vars &= !vars
    # @test !connected_path(_2_4, e3_4)
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

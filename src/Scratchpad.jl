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
    Vis.draw_dot(graph)
    subs, _ = compute_factorable_subgraphs(graph)

    _5_3 = subs[1]
    @assert (5, 3) == vertices(_5_3) #these are not tests. Put these here in case some future change to code causes order of subgraphs to change. Shouldn't happen but could.
    _2_4 = subs[2]
    @assert (2, 4) == vertices(_2_4)
    _3_5 = subs[4]
    @assert (3, 5) == vertices(_3_5)

    e5_3 = edges(graph, 5, 3)[1]

    pedges = collect(path(_5_3, e5_3))
    @assert length(pedges) == 1
    @assert e5_3 in pedges

    e3_4 = edges(graph, 3, 4)[1]
    e5_4 = edges(graph, 5, 4)[1]

    pedges = collect(path(_5_3, e3_4))
    @assert length(pedges) == 2
    @assert all(in.((e3_4, e5_4), Ref(pedges)))

    e2_3 = edges(graph, 2, 3)[1]
    e2_4 = edges(graph, 2, 4)[1]

    pedges = collect(path(_2_4, e3_4))
    @assert length(pedges) == 2
    @assert all(in.((e2_3, e3_4), Ref(pedges)))

    pedges = collect(path(2_4, e2_4))
    @assert length(pedges) == 1
    @assert e2_4 in pedges

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

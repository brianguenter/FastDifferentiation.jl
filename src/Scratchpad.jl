#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables x y
println("here")
    nv1 = Node(x)
    nv2 = Node(y)
    n3 = nv1 * nv2
    n4 = n3 * nv1

    n5 = n3 * n4

    graph = DerivativeGraph([n4, n5])
    # factor_subgraph!(graph, postdominator_subgraph(2, 4, 2, BitVector([0, 1]), BitVector([0, 1])))
    sub_heap, _ = compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)

    _5_3 = dominator_subgraph(graph, 5, 3, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _1_4 = postdominator_subgraph(graph, 1, 4, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _3_5 = postdominator_subgraph(graph, 3, 5, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _4_1 = dominator_subgraph(graph, 4, 1, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _5_1 = dominator_subgraph(graph, 5, 1, Bool[0, 1], Bool[0, 1], Bool[1, 0])
    _1_5 = postdominator_subgraph(graph, 1, 5, Bool[1, 0], Bool[0, 1], Bool[1, 0])

    correctly_ordered_subs = (_5_3, _1_4, _3_5, _4_1, _5_1, _1_5) #order of last two could switch and still be correct but all others should be in exactly this order.

    tmp = zip(correctly_ordered_subs[1:4], subs[1:4])
    for (correct, computed) in tmp
        @assert value_equal(correct, computed)
    end
    #last two
    @assert (value_equal(_5_1, subs[5]) && value_equal(_1_5, subs[6])) || (value_equal(_1_5, subs[5]) && value_equal(5_1, subs[6]))


end
export test


function bad_case()
    @variables x
    nx = Node(x)

    n2 = cos(nx)
    n3 = sin(nx)
    n4 = n2 * n3
    n5 = n2 * n4
    n6 = n2 - n3
    n7 = n5 * n6

    gr = DerivativeGraph(n7)
    println(map(x -> vertices(x), compute_factorable_subgraphs(gr)[1]))
    #if swap (1,4) and (7,3) in factorable subgraph list then get the error.
    #(1,4) comes first because its node difference is smaller than (7,3), which comes second. 
    #in this case can swap the order without introducing an ordering error. But the
    #factorization of (7,3) first creates a new factorable subgraph where factoring (1,4)
    #first didn't. 
    println(values(gr.postorder_number))
    Vis.draw_dot(gr, graph_label="before factoring")
    factor!(gr)
end
export bad_case

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

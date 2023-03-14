#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using Profile
using PProf
using .TestCases

function test()
    _, graph, four_2_subgraph, one_3_subgraph = simple_dominator_graph()



    subgraphs, subs_dict = compute_factorable_subgraphs(graph)
    @assert length(subgraphs) == 4
    four_one = subgraphs[findfirst(x -> x.subgraph == (4, 1), subgraphs)]
    one_3 = subgraphs[findfirst(x -> x.subgraph == (1, 3), subgraphs)]
    four_2 = subgraphs[findfirst(x -> x.subgraph == (4, 2), subgraphs)]
    one_4 = subgraphs[findfirst(x -> x.subgraph == (1, 4), subgraphs)]

    @assert factor_order(one_3, four_one) == true
    @assert factor_order(four_2, four_one) == true
    @assert factor_order(one_3, four_2) == false
    @assert factor_order(four_2, one_3) == false
    @assert factor_order(four_one, one_3) == false
    @assert factor_order(four_one, four_2) == false


    @assert factor_order(one_3, one_4) == true
    @assert factor_order(four_2, one_4) == true
    @assert factor_order(one_3, four_2) == false
    @assert factor_order(four_2, one_3) == false
    @assert factor_order(one_4, one_3) == false
    @assert factor_order(one_4, four_2) == false #one_3 should be sorted before one_4

    equal_subgraphs(x, y) = dominating_node(x) == dominating_node(y) && dominated_node(x) == dominated_node(y) && times_used(x) == times_used(y) && dominance_mask(x) == dominance_mask(y)


    # doms = dominator_subgraph.((
    #     (graph, 4, 2, BitVector([1, 0]), BitVector([1, 0])BitVector([1])),
    #     (graph, 4, 1, BitVector([1, 1]), BitVector([1, 0])BitVector([1]))))
    # pdoms = postdominator_subgraph.((
    #     (graph, 1, 3, BitVector([1, 1]), BitVector([1]), BitVector([1, 1])),
    #     (graph, 1, 4, BitVector([1, 1]), BitVector([1]), BitVector([1, 0]))))
    # subs2 = collect((pdoms..., doms...))


    index_1_4 = findfirst(x -> equal_subgraphs(x, one_4), subgraphs)
    index_4_2 = findfirst(x -> equal_subgraphs(x, four_2), subgraphs)
    index_1_3 = findfirst(x -> equal_subgraphs(x, one_3), subgraphs)

    @assert index_1_3 < index_1_4

end
export test

function profile()
    graph, qx, qy, qz = to_graph(25)
    Profile.clear()
    @profile symbolic_jacobian!(graph, [Node(qx), Node(qy), Node(qz)])
    pprof()
end
export profile

#changed
#change
export test

#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables x, y

    nx1 = Node(x)
    ny2 = Node(y)
    nxy3 = nx1 * ny2
    r2_4 = nx1 * nxy3
    r1_5 = r2_4 * nxy3

    gnodes = (nx1, ny2, nxy3, r2_4, r1_5)

    graph = DerivativeGraph([r1_5, r2_4])

    #first verify all nodes have the postorder numbers we expect
    for (i, nd) in pairs(gnodes)
        @assert node(graph, i) == nd
    end

    subs, subs_dict = compute_factorable_subgraphs(graph)
    sub_5_3 = first(filter(x -> x.subgraph == (5, 3), subs))
    sub_3_5 = first((filter(x -> x.subgraph == (3, 5), subs)))

    rmask = dominance_mask(sub_5_3)
    V = reachable_variables(sub_5_3)
    up_constraint = PathConstraint(dominating_node(sub_5_3), graph, true, rmask, V)
    dominating = sub_5_3.subgraph[1]
    dominated = sub_5_3.subgraph[2]
    path_edges1 = [edges(graph, 4, 3)[1], edges(graph, 5, 4)[1]]
    path_edges2 = [edges(graph, 5, 3)[1]]

    #test for dominator subgraph (5,3)
    up_edges = relation_edges!(up_constraint, dominated)

    temp_edges = PathEdge{Int64}[]

    edges_on_path!(up_constraint, dominating, true, up_edges[1], temp_edges)
    @assert all(x -> x[1] == x[2], zip(path_edges1, temp_edges))
    edges_on_path!(up_constraint, dominating, true, up_edges[2], temp_edges)
    @assert all(x -> x[1] == x[2], zip(path_edges2, temp_edges))

    #test for postdominator subgraph (3,5)
    down_constraint = PathConstraint(dominating_node(sub_3_5), graph, false, reachable_roots(sub_3_5), dominance_mask(sub_3_5))
    down_edges = relation_edges!(down_constraint, dominating)
    edges_on_path!(down_constraint, dominated, false, down_edges[1], temp_edges)
    @assert all(x -> x[1] == x[2], zip(reverse(path_edges1), temp_edges))
    edges_on_path!(down_constraint, dominated, false, down_edges[2], temp_edges)
    @assert all(x -> x[1] == x[2], zip(path_edges2, temp_edges))


    path_edges1 = [edges(graph, 4, 1)[1]]
    path_edges2 = [edges(graph, 3, 1)[1], edges(graph, 4, 3)[1]]
    sub_4_1 = first(filter(x -> x.subgraph == (4, 1), subs))
    sub_1_4 = first(filter(x -> x.subgraph == (1, 4), subs))
    dominating = sub_4_1.subgraph[1]
    dominated = sub_4_1.subgraph[2]

    #test for dominator subgraph (4,1)
    up_constraint = PathConstraint(dominating_node(sub_4_1), graph, true, dominance_mask(sub_4_1), reachable_variables(sub_4_1))
    up_edges = relation_edges!(up_constraint, dominated)
    edges_on_path!(up_constraint, dominating, true, up_edges[1], temp_edges)
    @assert all(x -> x[1] == x[2], zip(path_edges1, temp_edges))
    edges_on_path!(up_constraint, dominating, true, up_edges[2], temp_edges)
    @assert all(x -> x[1] == x[2], zip(path_edges2, temp_edges))

    dominating = sub_1_4.subgraph[1]
    dominated = sub_1_4.subgraph[2]
    #test for postdominator subgraph (1,4)
    down_constraint = PathConstraint(dominating_node(sub_1_4), graph, false, reachable_roots(sub_1_4), dominance_mask(sub_1_4))
    down_edges = relation_edges!(down_constraint, dominated)
    edges_on_path!(down_constraint, dominating, false, down_edges[1], temp_edges)
    @assert all(x -> x[1] == x[2], zip(path_edges1, temp_edges))
    edges_on_path!(down_constraint, dominating, false, down_edges[2], temp_edges)
    @assert all(x -> x[1] == x[2], zip(reverse(path_edges2), temp_edges))
end
export test

#changed
#change
export test

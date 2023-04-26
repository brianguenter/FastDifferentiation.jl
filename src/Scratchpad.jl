#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    dgraph = DerivativeGraph([complex_dominator_dag()])

    _1_4_sub_ref = Set(map(x -> x[1], edges.(Ref(dgraph), ((4, 3), (4, 2), (2, 1), (3, 1)))))

    _8_4_sub_ref = Set(map(x -> x[1], edges.(Ref(dgraph), ((8, 7), (8, 5), (5, 4), (7, 4)))))

    subs = extract_all!(compute_factorable_subgraphs(dgraph))
    _1_4_sub = subs[findfirst(x -> vertices(x) == (1, 4), subs)]
    _1_7_sub = subs[findfirst(x -> vertices(x) == (1, 7), subs)]
    _8_4_sub = subs[findfirst(x -> vertices(x) == (8, 4), subs)]
    _8_1_sub = subs[findfirst(x -> vertices(x) == (8, 1), subs)]
    _1_8_sub = subs[findfirst(x -> vertices(x) == (1, 8), subs)]

    @assert issetequal(_1_4_sub_ref, subgraph_edges(_1_4_sub))
    factor_subgraph!(_1_4_sub)
    _1_7_sub_ref = Set(map(x -> x[1], edges.(Ref(dgraph), ((4, 1), (3, 1), (7, 4), (7, 6), (6, 3)))))


    @assert issetequal(_1_7_sub_ref, subgraph_edges(_1_7_sub))
    @assert issetequal(_8_4_sub_ref, subgraph_edges(_8_4_sub))
    factor_subgraph!(_8_4_sub)
    _8_1_sub_ref = Set(map(x -> x[1], edges.(Ref(dgraph), ((8, 7), (8, 4), (4, 1), (3, 1), (6, 3), (7, 6)))))
    @assert issetequal(_8_1_sub_ref, subgraph_edges(_8_1_sub))
    @assert issetequal(_8_1_sub_ref, subgraph_edges(_1_8_sub))

end
export test


function test2()
    @variables x y

    nx = Node(x)
    func = nx * nx

    gr = DerivativeGraph([func])
    subs_heap = compute_factorable_subgraphs(gr)
    subs = extract_all!(subs_heap)
    test_sub = subs[1]
    etmp = parent_edges(gr, dominated_node(test_sub))
    rroots = reachable_roots(etmp[1])
    rroots .= rroots .& .!rroots

    @assert !connected_path(test_sub, etmp[1])
    @assert connected_path(test_sub, etmp[2])
end
export test2

function make_relations(size)
    relations = Vector{Vector{Int64}}(undef, size)
    for i in eachindex(relations)
        if i == size
            relations[i] = Int64[]
        elseif i == size - 1
            relations[i] = [size]
        elseif i == size - 2
            relations[i] = [size - 1, size]
        else
            relations[i] = [i + 1, i + 2]
        end
    end
    return relations
end
export make_relations

function time_relations()
    gr = simple_dominator_dgraph()[1]
    dp = DomPathConstraint(gr, true, 1)
    return gr, dp
end
export time_relations

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

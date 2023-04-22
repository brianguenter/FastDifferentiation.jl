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
    subs_heap, _ = compute_factorable_subgraphs(gr)
end
export test

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

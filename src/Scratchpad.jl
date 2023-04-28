#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    graph, subs = simple_factorable_subgraphs()

    all_edges = collect(unique_edges(graph))

    _4_2 = all_edges[findfirst(x -> vertices(x) == (4, 2), all_edges)]
    _4_3 = all_edges[findfirst(x -> vertices(x) == (4, 3), all_edges)]
    _3_2 = all_edges[findfirst(x -> vertices(x) == (3, 2), all_edges)]
    _2_1 = all_edges[findfirst(x -> vertices(x) == (2, 1), all_edges)]
    _3_1 = all_edges[findfirst(x -> vertices(x) == (3, 1), all_edges)]

    ed, nod = deconstruct_subgraph(subs[1]) #can only deconstruct these two subgraphs because the larger ones need to be factored first.
    @assert issetequal([4, 2, 3], nod)
    @assert issetequal((_4_2, _4_3, _3_2), ed)

    ed, nod = deconstruct_subgraph(subs[2])
    @assert issetequal((_3_2, _3_1, _2_1), ed)
    @assert issetequal([1, 2, 3], nod)

    factor_subgraph!(subs[1]) #now can test larger subgraphs

    #new edges created during factorization so get them again
    all_edges = collect(unique_edges(graph))

    _4_2 = all_edges[findfirst(x -> vertices(x) == (4, 2), all_edges)]
    _4_3 = all_edges[findfirst(x -> vertices(x) == (4, 3), all_edges)]
    _2_1 = all_edges[findfirst(x -> vertices(x) == (2, 1), all_edges)]
    _3_1 = all_edges[findfirst(x -> vertices(x) == (3, 1), all_edges)]

    ed, nod = deconstruct_subgraph(subs[3])
    println(ed)
    sub_4_1 = (_4_3, _4_2, _3_1, _2_1)
    @assert issetequal(sub_4_1, ed)
    @assert issetequal([1, 2, 3, 4], nod)
    ed, nod = deconstruct_subgraph(subs[4])
    @assert issetequal(sub_4_1, ed)
    @assert issetequal([1, 2, 3, 4], nod)
end
export make_subdata

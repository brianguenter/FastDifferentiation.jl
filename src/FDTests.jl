#Test functions need access to non-exported names from FastDifferentiation. 
module FDInternals
using FastDifferentiation: compute_factorable_subgraphs, edges, vertices, isa_connected_path, add_non_dom_edges!, edges, PathEdge, dominated_node, partial_value, compute_paths_to_roots, num_vertices, unique_edges, parent_edges, edge_path, all_nodes, dominator_subgraph, postdominator_subgraph, dominating_node, times_used, reachable_dominance, postorder_number, variable_index_to_postorder_number, root_index_to_postorder_number, parents, children, DomPathConstraint, edge_exists, FactorableSubgraph, each_vertex, node_edges, factor_order, subgraph_edges, node, subset, factor_subgraph!, factor!, deconstruct_subgraph, forward_edges, evaluate_subgraph, make_factored_edge, path_sort_order, multiply_sequence, root, reachable_roots, bott_vertex, top_vertex, reachable_variables, compute_paths_to_variables, compute_edge_paths!, relation_node_indices, value_equal, is_tree, is_leaf, is_variable, is_constant, set_diff, value, bit_equal, _symbolic_jacobian, _symbolic_jacobian!, _sparse_symbolic_jacobian!, roots, variables, domain_dimension, codomain_dimension, DerivativeGraph, Node

export compute_factorable_subgraphs, edges, vertices, isa_connected_path, add_non_dom_edges!, edges, PathEdge, dominated_node, partial_value, compute_paths_to_roots, num_vertices, unique_edges, parent_edges, edge_path, all_nodes, dominator_subgraph, postdominator_subgraph, dominating_node, times_used, reachable_dominance, postorder_number, variable_index_to_postorder_number, root_index_to_postorder_number, parents, children, DomPathConstraint, edge_exists, FactorableSubgraph, each_vertex, node_edges, factor_order, subgraph_edges, node, subset, factor_subgraph!, factor!, deconstruct_subgraph, forward_edges, evaluate_subgraph, make_factored_edge, path_sort_order, multiply_sequence, root, reachable_roots, bott_vertex, top_vertex, reachable_variables, compute_paths_to_variables, compute_edge_paths!, relation_node_indices, value_equal, is_tree, is_leaf, is_variable, is_constant, set_diff, value, bit_equal, _symbolic_jacobian, _symbolic_jacobian!, _sparse_symbolic_jacobian!, roots, variables, domain_dimension, codomain_dimension, DerivativeGraph, Node
end #module


module FDTests
using DiffRules
using TermInterface
import SymbolicUtils
using StaticArrays
using Memoize
using DataStructures

using ..FastDifferentiation
using ..FDInternals

const FD = FastDifferentiation
export FD

using TestItems

include("../TestPrograms/TestCode.jl")

"""If `compute_dominators` is `true` then computes `idoms` tables for graph, otherwise computes `pidoms` table`"""
function compute_dominance_tables(graph::DerivativeGraph{T}, compute_dominators::Bool) where {T<:Integer}
    if compute_dominators
        start_vertices = root_index_to_postorder_number(graph)
    else
        start_vertices = variable_index_to_postorder_number(graph)
    end

    doms = Dict{T,T}[]   #create one idom table for each root

    for (start_index, node_postorder_number) in pairs(start_vertices)
        push!(doms, FastDifferentiation.compute_dom_table(graph, compute_dominators, start_index, node_postorder_number))
    end
    return doms
end
export compute_dominance_tables

""" Utility function for working with symbolic expressions as Symbolics.jl defines them."""
function number_of_operations(symbolic_expr)
    if SymbolicUtils.istree(symbolic_expr) && operation(symbolic_expr) ∈ (+, *, -)
        return 1 + sum(number_of_operations.(arguments(symbolic_expr)))
    else
        return 0
    end
end

function simple_dag(cache::Union{IdDict,Nothing}=IdDict())
    @variables zz y
    return expr_to_dag(zz^2 + y * (zz^2), cache), zz, y
end
export simple_dag

function simple_numbered_dag()
    @variables zz
    return expr_to_dag(zz * cos(zz^2))
end
export simple_numbered_dag

function dominators_dag()
    @variables zz
    return expr_to_dag(zz * (cos(zz) + sin(zz))), zz
end
export dominators_dag

#Each line in the postorder listing has two spaces at the end which causes a line break. Don't delete the trailing spaces or the formatting will get messed up.
"""
Every line in this comment is terminated by 2 spaces to make markdown format properly.

creates dag with this postorder:  
1    x  
2    cos (x)  
3    sin (x)  
4    k1 = (cos(x) * sin(x))  
5    sin (k1)  
6    k2 = exp(sin(x))  
7    (k1 * k2)  
8    (sin(k1) + (k1 * k2))  

with this edge table:  
nodenum  
1    [(2, 1), (3, 1)]  
2    [(2, 1), (4, 2)]  
3    [(3, 1), (4, 3), (6, 3)]  
4    [(4, 2), (4, 3), (5, 4), (7, 4)]  
5    [(5, 4), (8, 5)]  
6    [(6, 3), (7, 6)]  
7    [(7, 4), (7, 6), (8, 7)]  
8    [(8, 5), (8, 7)]  

with this idoms/pidoms table:  
nodenum   idoms    pidoms  
1         8        1  
2         4        1  
3         8        1  
4         8        1  
5         8        4  
6         7        3  
7         8        1  
8         8        1  

with these factorable subgraphs  
(8,4), (8,3), (8,1), (1,4), (1,7), (1,8)
"""
function complex_dominator_dag()
    @variables x

    sinx = FastDifferentiation.Node(sin, MVector(x))
    cosx = FastDifferentiation.Node(cos, MVector(x))
    A = FastDifferentiation.Node(*, MVector(cosx, sinx))
    sinA = FastDifferentiation.Node(sin, MVector(A))
    expsin = FastDifferentiation.Node(*, MVector(A, FastDifferentiation.Node(exp, MVector(sinx))))
    plus = FastDifferentiation.Node(+, MVector(sinA, expsin))
    return plus
end
export complex_dominator_dag

complex_dominator_graph() = DerivativeGraph(complex_dominator_dag())
export complex_dominator_graph

function R2_R2_function()
    @variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4

    return DerivativeGraph([n5, n4])
end
export R2_R2_function

function simple_dominator_graph()
    @variables x

    ncos = Node(cos, x)
    nplus = Node(+, ncos, x)
    ntimes = Node(*, ncos, nplus)
    four_2_subgraph = Node(+, nplus, ncos)
    one_3_subgraph = Node(+, Node(*, Node(-1), Node(sin, x)), Node(1))
    return x, DerivativeGraph(ntimes), four_2_subgraph, one_3_subgraph
end
export simple_dominator_graph

"""returns 4 factorable subgraphs in this order: (4,2),(1,3),(1,4),(4,1)"""
function simple_factorable_subgraphs()
    _, graph, _, _ = simple_dominator_graph()
    temp = extract_all!(compute_factorable_subgraphs(graph))
    return graph, [
        temp[findfirst(x -> FastDifferentiation.vertices(x) == (4, 2), temp)],
        temp[findfirst(x -> FastDifferentiation.vertices(x) == (1, 3), temp)],
        temp[findfirst(x -> FastDifferentiation.vertices(x) == (1, 4), temp)],
        temp[findfirst(x -> FastDifferentiation.vertices(x) == (4, 1), temp)]
    ]
end
export simple_factorable_subgraphs
@testitem "isa_connected_path 1" begin # case when path is one edge long
    using DataStructures
    using FastDifferentiation.FDInternals

    @variables x y

    func = x * x

    gr = DerivativeGraph([func])
    subs_heap = compute_factorable_subgraphs(gr)
    subs = extract_all!(subs_heap)
    test_sub = subs[1]

    etmp = parent_edges(gr, dominated_node(test_sub))
    rroots = reachable_roots(etmp[1])
    rroots .= rroots .& .!rroots

    @test !isa_connected_path(test_sub, etmp[1])
    @test isa_connected_path(test_sub, etmp[2])
end

@testitem "isa_connected_path 2" begin #cases when path is longer than one edge and various edges have either roots or variables reset.

    using DataStructures
    using FastDifferentiation.FDInternals
    @variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4


    graph = DerivativeGraph([n4, n5])
    subs_heap = compute_factorable_subgraphs(graph)
    subs = extract_all!(subs_heap)
    println(subs)
    _5_3_index = findfirst(x -> vertices(x) == (5, 3), subs)
    _5_3 = subs[_5_3_index]

    _2_4_index = findfirst(x -> vertices(x) == (2, 4), subs)
    _2_4 = subs[_2_4_index]

    _3_5_index = findfirst(x -> vertices(x) == (3, 5), subs)
    _3_5 = subs[_3_5_index]

    etmp = edges(graph, 3, 5)[1]
    @test isa_connected_path(_5_3, etmp)


    etmp = edges(graph, 3, 4)[1]
    @test isa_connected_path(_5_3, etmp)
    rts = reachable_roots(etmp)
    rts[2] = 0

    @test !isa_connected_path(_5_3, etmp)
    #reset path
    rts[2] = 1

    e2_4 = edges(graph, 2, 4)[1]
    @test isa_connected_path(_2_4, e2_4)
    e2_3 = edges(graph, 2, 3)[1]
    @test isa_connected_path(_2_4, e2_3)
    e3_4 = edges(graph, 3, 4)[1]
    vars = reachable_variables(e3_4)
    @. vars &= !vars
    @test !isa_connected_path(_2_4, e3_4)
end

@testitem "add_non_dom_edges" begin

    using DataStructures
    using FastDifferentiation.FDInternals

    #utility function to make it easier to create edges and test them against edges generated during graph operations.
    function edge_fields_equal(edge1, edge2)
        return edge1.top_vertex == edge2.top_vertex &&
               edge1.bott_vertex == edge2.bott_vertex &&
               edge1.edge_value == edge2.edge_value &&
               edge1.reachable_variables == edge2.reachable_variables &&
               edge1.reachable_roots == edge2.reachable_roots
    end

    @variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4

    graph = DerivativeGraph([n4, n5])
    subs_heap = compute_factorable_subgraphs(graph)
    subs = extract_all!(subs_heap)
    _5_3 = subs[1]
    @test (5, 3) == vertices(_5_3)

    add_non_dom_edges!(_5_3)
    #single edge 3,4 should be split into two: ([r1,r2],[v1,v2]) -> ([r1],[v1,v2]),([r2],[v1,v2])
    edges3_4 = edges(graph, 4, 3)
    @test length(edges3_4) == 2
    test_edge = PathEdge(4, 3, y, BitVector([1, 1]), BitVector([0, 1]))
    @test count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1
    test_edge = (PathEdge(4, 3, y, BitVector([1, 1]), BitVector([1, 0])))
    @test count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1

    graph = DerivativeGraph([n4, n5])
    sub_heap = compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)
    _2_4 = subs[2]
    @test (2, 4) == vertices(_2_4)

    add_non_dom_edges!(_2_4)
    #single edge 3,4 should be split in two: ([r1,r2],[v1,v2])->([r1,r2],[v1]),([r1,r2],[v2])
    edges3_4 = edges(graph, 4, 3)
    @test length(edges3_4) == 2
    test_edge = PathEdge(4, 3, y, BitVector([1, 0]), BitVector([1, 1]))
    @test count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1
    test_edge = (PathEdge(4, 3, y, BitVector([0, 1]), BitVector([1, 1])))
    @test count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1
end

@testitem "iteration" begin

    using DataStructures
    using FastDifferentiation.FDInternals

    @variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4


    graph = DerivativeGraph([n4, n5])
    subs_heap = compute_factorable_subgraphs(graph)

    subs = extract_all!(subs_heap)

    _5_3_index = findfirst(x -> vertices(x) == (5, 3), subs)
    _5_3 = subs[_5_3_index]

    _2_4_index = findfirst(x -> vertices(x) == (2, 4), subs)
    _2_4 = subs[_2_4_index]

    _3_5_index = findfirst(x -> vertices(x) == (3, 5), subs)
    _3_5 = subs[_3_5_index]

    e5_3 = edges(graph, 5, 3)[1]

    pedges = collect(edge_path(_5_3, e5_3))
    @test length(pedges) == 1
    @test e5_3 in pedges

    e3_4 = edges(graph, 3, 4)[1]
    e5_4 = edges(graph, 5, 4)[1]

    pedges = collect(edge_path(_5_3, e3_4))
    @test length(pedges) == 2
    @test all(in.((e3_4, e5_4), Ref(pedges)))

    e2_3 = edges(graph, 2, 3)[1]
    e2_4 = edges(graph, 2, 4)[1]

    pedges = collect(edge_path(_2_4, e3_4))
    @test length(pedges) == 2
    @test all(in.((e2_3, e3_4), Ref(pedges)))

    pedges = collect(edge_path(_2_4, e2_4))
    @test length(pedges) == 1
    @test e2_4 in pedges
end



@testitem "is_tree" begin
    using FastDifferentiation.FDInternals

    @variables x

    x = Node(x)
    z = Node(0)
    tm = x * x

    @test is_tree(x) == false
    @test is_leaf(x) == true
    @test is_variable(x) == true
    @test is_constant(x) == false

    @test is_tree(z) == false
    @test is_leaf(z) == true
    @test is_constant(z) == true
    @test is_variable(z) == false

    @test is_leaf(tm) == false
    @test is_variable(tm) == false
    @test is_tree(tm) == true
    @test is_constant(tm) == false
end

@testitem "derivative" begin
    using FastDifferentiation.FDInternals
    @variables x y


    a = x * y
    @test derivative(a, Val(1)) == y
    @test derivative(a, Val(2)) == x
    @test derivative(x) == Node(1)
    @test derivative(Node(1)) == Node(0)
end

@testitem "compute_factorable_subgraphs test order" begin

    using DataStructures
    using FastDifferentiation.FDInternals

    @variables x y

    nv1 = Node(x)
    nv2 = Node(y)
    n3 = nv1 * nv2
    n4 = n3 * nv1

    n5 = n3 * n4

    graph = DerivativeGraph([n4, n5])
    # factor_subgraph!(graph, postdominator_subgraph(2, 4, 2, BitVector([0, 1]), BitVector([0, 1])))
    sub_heap = compute_factorable_subgraphs(graph)
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
        @test value_equal(correct, computed)
    end
    #last two
    @test (value_equal(_5_1, subs[5]) && value_equal(_1_5, subs[6])) || (value_equal(_1_5, subs[5]) && value_equal(5_1, subs[6]))
end

@testitem "compute_factorable_subgraphs" begin
    using FastDifferentiation.FDTests
    using DataStructures
    using FastDifferentiation.FDInternals

    dgraph = DerivativeGraph(complex_dominator_dag())

    sub_heap = compute_factorable_subgraphs(dgraph)
    subs = extract_all!(sub_heap)


    equal_subgraphs(x, y) = dominating_node(x) == dominating_node(y) && dominated_node(x) == dominated_node(y) && times_used(x) == times_used(y) && reachable_dominance(x) == reachable_dominance(y)


    index_1_4 = findfirst(x -> equal_subgraphs(x, postdominator_subgraph(dgraph, 1, 4, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    index_1_7 = findfirst(x -> equal_subgraphs(x, postdominator_subgraph(dgraph, 1, 7, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    @test index_1_4 < index_1_7

    index_8_4 = findfirst(x -> equal_subgraphs(x, dominator_subgraph(dgraph, 8, 4, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    index_8_3 = findfirst(x -> equal_subgraphs(x, dominator_subgraph(dgraph, 8, 3, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    @test index_8_4 < index_8_3

    index_8_1d = findfirst(x -> equal_subgraphs(x, dominator_subgraph(dgraph, 8, 1, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    @test index_8_4 < index_8_1d
    @test index_8_3 < index_8_1d
    @test index_1_7 < index_8_1d

    index_8_1p = findfirst(x -> equal_subgraphs(x, dominator_subgraph(dgraph, 8, 1, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    @test index_8_4 < index_8_1p
    @test index_8_3 < index_8_1p
    @test index_1_7 < index_8_1p
end

@testitem "edges" begin

    using FastDifferentiation.FDInternals

    @variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4

    graph = DerivativeGraph([n5, n4])

    function test_edge_access(graph, correct_num_edges, vert1, vert2)
        edge_group1 = edges(graph, vert1, vert2)
        edge_group2 = edges(graph, vert2, vert1)
        @test length(edge_group1) == correct_num_edges
        @test length(edge_group1) == length(edge_group2)
        for edge_pair in zip(edge_group1, edge_group2)
            @test edge_pair[1] == edge_pair[2]
        end

        for edge in edge_group1
            @test top_vertex(edge) == max(vert1, vert2)
            @test bott_vertex(edge) == min(vert1, vert2)
        end
    end

    edge_verts = ((1, 3), (3, 2), (2, 4), (5, 4), (3, 5), (4, 3))
    num_edges = (1, 1, 1, 1, 1, 1)

    for (edge, num) in zip(edge_verts, num_edges)
        test_edge_access(graph, num, edge[1], edge[2])
    end
end

@testitem "DerivativeGraph constructor" begin
    using FastDifferentiation.FDInternals

    @variables x y


    cosx = Node(cos, x)
    sinx = Node(sin, x)
    partial_cosx = Node(-, sinx)
    siny = Node(sin, y)
    partial_siny = Node(cos, y)
    ctimess = cosx * siny
    partial_times = [siny, cosx]
    cpluss = cosx + siny
    partial_plus = [Node(1), Node(1)]
    roots = [ctimess, cpluss]
    grnodes = [x, y, cosx, siny, cpluss, ctimess]

    correct_postorder_numbers = Dict((x => 1, cosx => 2, y => 3, siny => 4, ctimess => 5, cpluss => 6))

    graph = DerivativeGraph(roots)

    @test all([correct_postorder_numbers[node] == postorder_number(graph, node) for node in grnodes])

    correct_partials = Dict((cosx => [partial_cosx], siny => [partial_siny], ctimess => partial_times, cpluss => partial_plus))
    #TODO need to use finite differences instead of Symbolics. 
    # for (node, partials) in pairs(correct_partials)
    #     for (i, one_partial) in pairs(partials)
    #         f1 = to_symbolics(partial_value(graph, node, i))
    #         f2 = to_symbolics(one_partial)

    #         for test_point in BigFloat(-1):BigFloat(0.01):BigFloat(1) #graphs might have equivalent but different forms so evaluate at many points at high precision to verify derivatives are the same.
    #             v1 = Symbolics.value(Symbolics.substitute(f1, Dict((x => test_point), (y => test_point))))
    #             v2 = Symbolics.value(Symbolics.substitute(f2, Dict((x => test_point), (y => test_point))))
    #             @test isapprox(v1, v2, atol=1e-50)
    #         end
    #     end
    # end
end

@testitem "DerivativeGraph pathmasks" begin

    using FastDifferentiation.FDInternals

    @variables x y

    xy = Node(*, x, y) #postorder # 4
    n5 = Node(5) #postorder # 1
    f1 = Node(*, n5, xy) #postorder 5
    n3 = Node(3) #postorder # 6
    f2 = Node(*, xy, n3) #postorder # 7
    roots = [f1, f2]
    graph = DerivativeGraph(roots)

    #x,y,xy shared by both roots. n5,f1 only on f1 path, n3,f2 only on f2 path
    correct_roots_pathmasks = [
        BitArray([1, 0]),
        BitArray([1, 1]),
        BitArray([1, 1]),
        BitArray([1, 1]),
        BitArray([1, 0]),
        BitArray([0, 1]),
        BitArray([0, 1])]

    correct_variable_pathmasks = [
        BitArray([0, 0]),
        BitArray([1, 0]),
        BitArray([0, 1]),
        BitArray([1, 1]),
        BitArray([1, 1]),
        BitArray([0, 0]),
        BitArray([1, 1])
    ]

    variable_path_masks = compute_paths_to_variables(num_vertices(graph), edges(graph), variable_index_to_postorder_number(graph))
    @test variable_path_masks == correct_variable_pathmasks
    parent_path_masks = compute_paths_to_roots(num_vertices(graph), edges(graph), root_index_to_postorder_number(graph))
    @test parent_path_masks == correct_roots_pathmasks
end

@testitem "ConstrainedPathIterator" begin

    using FastDifferentiation.FDInternals

    @variables x y

    #ℝ²->ℝ² function (f1,f2) = (5*(x*y),(x*y)*3)

    xy = Node(*, x, y) #postorder # 4
    n5 = Node(5) #postorder # 1
    f1 = Node(*, n5, xy) #postorder 5
    n3 = Node(3) #postorder # 6
    f2 = Node(*, xy, n3) #postorder # 7
    roots = [f1, f2]
    graph = DerivativeGraph(roots)
    root_masks = compute_paths_to_roots(num_vertices(graph), edges(graph), root_index_to_postorder_number(graph))
    variable_masks = compute_paths_to_variables(num_vertices(graph), edges(graph), variable_index_to_postorder_number(graph))

    piterator = DomPathConstraint(graph, true, 1)

    parents = Int64[]
    correct_parents = (postorder_number(graph, xy))
    for parent in relation_node_indices(piterator, postorder_number(graph, x))
        push!(parents, parent)
    end

    @test length(parents) == 1
    @test parents[1] == postorder_number(graph, xy)

    piterator = DomPathConstraint(graph, true, 1)
    parents = Int64[]
    pnum = postorder_number(graph, xy)
    for parent in relation_node_indices(piterator, pnum)
        push!(parents, parent)
    end
    @test length(parents) == 1
    @test postorder_number(graph, f1) == parents[1]

    piterator = DomPathConstraint(graph, true, 2)

    parents = Int64[]
    for parent in relation_node_indices(piterator, postorder_number(graph, xy))
        push!(parents, parent)
    end

    @test length(parents) == 1
    @test postorder_number(graph, f2) == parents[1]


    viterator = DomPathConstraint(graph, false, 1)
    children = Int64[]
    for child in relation_node_indices(viterator, postorder_number(graph, xy))
        push!(children, child)
    end
    @test length(children) == 1
    @test children[1] == postorder_number(graph, x)

    viterator = DomPathConstraint(graph, false, 2)
    children = Int64[]
    for child in relation_node_indices(viterator, postorder_number(graph, xy))
        push!(children, child)
    end
    @test length(children) == 1
    @test children[1] == 3
end

@testitem "edge_exists" begin

    using FastDifferentiation.FDInternals

    @variables x y

    xy = Node(*, x, y) #postorder # 4
    n5 = Node(5) #postorder # 1
    f1 = Node(*, n5, xy) #postorder 5
    n3 = Node(3) #postorder # 6
    f2 = Node(*, xy, n3) #postorder # 7
    roots = [f1, f2]
    graph = DerivativeGraph(roots)

    all_edges = unique_edges(graph)

    for edge in all_edges
        @test edge_exists(graph, edge)
    end

    @test !edge_exists(graph, PathEdge(1, 7, Node(0), 2, 2)) #this edge is not in the graph

end

@testitem "add_edge! for DerivativeGraph" begin

    using FastDifferentiation.FDInternals

    @variables x y

    xy = Node(*, x, y) #postorder # 4
    n5 = Node(5) #postorder # 1
    f1 = Node(*, n5, xy) #postorder 5
    n3 = Node(3) #postorder # 6
    f2 = Node(*, xy, n3) #postorder # 7
    roots = [f1, f2]
    variables = [x, y]
    graph = DerivativeGraph(roots)

    previous_edges = unique_edges(graph)
    new_edge = PathEdge(1, 7, Node(y), length(variables), length(roots))
    FastDifferentiation.add_edge!(graph, new_edge)

    #make sure existing edges are still in the graph.
    for edge in previous_edges
        @test edge_exists(graph, edge)
    end

    prnts = parents.(values(edges(graph)))
    numprnts = sum(length.(prnts))
    chldrn = children.(values(edges(graph)))
    numchldrn = sum(length.(chldrn))
    num_edges = (numprnts + numchldrn) / 2
    @test num_edges == 7 #ensure number of edges has increased by 1

    @test edge_exists(graph, new_edge) #and that there is only one new edge
end

@testitem "delete_edge! for DerivativeGraph" begin

    using FastDifferentiation.FDInternals

    @variables x y

    function reset_test(all_edges, graph, func::Function)
        for edge in all_edges
            tmp = func(edge)
            for i in eachindex(tmp)
                tmp[i] = 0
            end

            #can't delete edge till roots or variables are all unreachable

            FastDifferentiation.delete_edge!(graph, edge)
            @test !edge_exists(graph, edge) #make sure edge has been deleted from graph

            delete!(all_edges, edge) #now delete edge and see if all the other edges that are still supposed to be in the graph are still there
            for edge2 in all_edges
                @test edge_exists(graph, edge2) #other edges have not been deleted
            end
        end
    end

    xy = Node(*, x, y) #postorder # 4
    n5 = Node(5) #postorder # 1
    f1 = Node(*, n5, xy) #postorder 5
    n3 = Node(3) #postorder # 6
    f2 = Node(*, xy, n3) #postorder # 7
    roots = [f1, f2]
    graph = DerivativeGraph(roots)
    all_edges = unique_edges(graph)

    reset_test(all_edges, graph, reachable_roots)

    graph = DerivativeGraph(roots)
    all_edges = unique_edges(graph)

    reset_test(all_edges, graph, reachable_variables)
end

@testitem "compute_edge_paths" begin

    using FastDifferentiation.FDInternals


    @variables x y

    xy = Node(*, x, y) #postorder # 4
    n5 = Node(5) #postorder # 1
    f1 = Node(*, n5, xy) #postorder 5
    n3 = Node(3) #postorder # 6
    f2 = Node(*, xy, n3) #postorder # 7
    roots = [f1, f2]
    variables = [x, y]
    graph = DerivativeGraph(roots)
    compute_edge_paths!(num_vertices(graph), edges(graph), variable_index_to_postorder_number(graph), root_index_to_postorder_number(graph))

    correct_root_masks = Dict(
        (4, 2) => BitVector([1, 1]),
        (4, 3) => BitVector([1, 1]),
        (5, 1) => BitVector([1, 0]),
        (5, 4) => BitVector([1, 0]),
        (7, 4) => BitVector([0, 1]),
        (7, 6) => BitVector([0, 1])
    )

    correct_variable_masks = Dict(
        (4, 2) => BitVector([1, 0]),
        (4, 3) => BitVector([0, 1]),
        (5, 1) => BitVector([0, 0]),
        (5, 4) => BitVector([1, 1]),
        (7, 4) => BitVector([1, 1]),
        (7, 6) => BitVector([0, 0])
    )

    for index in each_vertex(graph)
        c_and_p = node_edges(graph, index)
        for edge in [parents(c_and_p); children(c_and_p)]
            @test edge.reachable_variables == correct_variable_masks[(top_vertex(edge), bott_vertex(edge))]
            @test edge.reachable_roots == correct_root_masks[(top_vertex(edge), bott_vertex(edge))]
        end
    end
end

@testitem "dominators DerivativeGraph" begin

    using FastDifferentiation.FDTests
    using FastDifferentiation.FDInternals


    @variables x y

    xy = Node(*, x, y) #postorder # 4
    n5 = Node(5) #postorder # 1
    f1 = Node(*, n5, xy) #postorder 5
    n3 = Node(3) #postorder # 6
    f2 = Node(*, xy, n3) #postorder # 7
    roots = [f1, f2]
    graph = DerivativeGraph(roots)
    idoms = compute_dominance_tables(graph, true)

    correct_dominators = [
        (1 => 5, 4 => 5, 2 => 4, 3 => 4, 5 => 5),
        (2 => 4, 3 => 4, 6 => 7, 4 => 7, 7 => 7)
    ]

    for (i, idom) in pairs(idoms)
        for elt in correct_dominators[i]
            @test elt[2] == idom[elt[1]]
        end
    end


    ncos = Node(cos, x) # 2
    nsin = Node(sin, y) # 4
    n5 = Node(*, ncos, nsin) #5
    n6 = Node(*, n5, y) #6
    n7 = Node(*, n5, n6) #7
    nexp = Node(exp, n6) # 8

    roots = [n7, nexp]
    graph = DerivativeGraph(roots)
    idoms = compute_dominance_tables(graph, true)

    correct_dominators = [
        (1 => 2, 2 => 5, 3 => 7, 4 => 5, 5 => 7, 6 => 7, 7 => 7),
        (1 => 2, 2 => 5, 3 => 6, 4 => 5, 5 => 6, 6 => 8, 8 => 8)
    ]

    for (i, idom) in pairs(idoms)
        for elt in correct_dominators[i]
            @test elt[2] == idom[elt[1]]
        end
    end

    pidoms = compute_dominance_tables(graph, false)

    correct_post_dominators = [
        (1 => 1, 2 => 1, 5 => 2, 6 => 5, 7 => 5, 8 => 6),
        (3 => 3, 4 => 3, 5 => 4, 6 => 3, 7 => 3, 8 => 6)
    ]

    for (i, pidom) in pairs(pidoms)
        for elt in correct_post_dominators[i]
            @test elt[2] == pidom[elt[1]]
        end
    end
end


@testitem "dom_subgraph && pdom_subgraph" begin


    using FastDifferentiation: dom_subgraph, pdom_subgraph
    using FastDifferentiation.FDTests
    using FastDifferentiation.FDInternals

    @variables x y

    ncos = Node(cos, x) #pos
    ntimes1 = Node(*, ncos, x)
    ntimes2 = Node(*, ncos, ntimes1)
    roots = [ntimes2, ntimes1]
    graph = DerivativeGraph(roots)
    idoms = compute_dominance_tables(graph, true)
    pidoms = compute_dominance_tables(graph, false)

    r1_doms = ((4, 1), (4, 2))
    v1_pdoms = ((1, 3), (1, 4))

    computed = (
        dom_subgraph(graph, 1, 1, idoms[1]),
        dom_subgraph(graph, 1, 2, idoms[1]))

    @test computed[1] == r1_doms[1]
    @test computed[2] == r1_doms[2]

    computed = (
        pdom_subgraph(graph, 1, 3, pidoms[1]),
        pdom_subgraph(graph, 1, 4, pidoms[1]))

    @test computed[1] == v1_pdoms[1]
    @test computed[2] == v1_pdoms[2]

    r2_dom = (3, 1)

    computed = dom_subgraph(graph, 2, 1, idoms[2])

    @test computed == r2_dom
end


@testitem "reachable" begin

    using FastDifferentiation.FDInternals

    @variables x y


    nx1 = Node(x)
    ny2 = Node(y)
    n3 = nx1 * ny2
    n4 = n3 * ny2
    n5 = n3 * n4

    graph = DerivativeGraph([n5, n4])
    two = BitVector([1, 1])
    one_zero = BitVector([1, 0])
    zero_one = BitVector([0, 1])

    for node in (1, 2, 3, 4)
        @test reachable_roots(graph, node) == two
    end
    @test reachable_roots(graph, 5) == one_zero

    @test reachable_variables(graph, 2) == zero_one
    for node in (3, 4, 5)
        @test reachable_variables(graph, node) == two
    end

    @test reachable_variables(graph, 1) == one_zero
end

@testitem "relation_edges" begin

end

@testitem "factor_order" begin
    using FastDifferentiation.FDTests
    using DataStructures
    using FastDifferentiation.FDInternals

    _, graph, four_2_subgraph, one_3_subgraph = simple_dominator_graph()



    sub_heap = compute_factorable_subgraphs(graph)
    subgraphs = extract_all!(sub_heap)
    @test length(subgraphs) == 4
    four_one = subgraphs[findfirst(x -> x.subgraph == (4, 1), subgraphs)]
    one_3 = subgraphs[findfirst(x -> x.subgraph == (1, 3), subgraphs)]
    four_2 = subgraphs[findfirst(x -> x.subgraph == (4, 2), subgraphs)]
    one_4 = subgraphs[findfirst(x -> x.subgraph == (1, 4), subgraphs)]

    @test factor_order(one_3, four_one) == true
    @test factor_order(four_2, four_one) == true
    @test factor_order(one_3, four_2) == false
    @test factor_order(four_2, one_3) == false
    @test factor_order(four_one, one_3) == false
    @test factor_order(four_one, four_2) == false


    @test factor_order(one_3, one_4) == true
    @test factor_order(four_2, one_4) == true
    @test factor_order(one_3, four_2) == false
    @test factor_order(four_2, one_3) == false
    @test factor_order(one_4, one_3) == false
    @test factor_order(one_4, four_2) == false #one_3 should be sorted before one_4

    equal_subgraphs(x, y) = dominating_node(x) == dominating_node(y) && dominated_node(x) == dominated_node(y) && times_used(x) == times_used(y) && reachable_dominance(x) == reachable_dominance(y)


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

    @test index_1_3 < index_1_4
end

@testitem "subgraph_edges" begin
    using FastDifferentiation.FDTests
    using DataStructures
    using FastDifferentiation.FDInternals

    dgraph = DerivativeGraph([complex_dominator_dag()])

    _1_4_sub_ref = Set(map(x -> x[1], edges.(Ref(dgraph), ((4, 3), (4, 2), (2, 1), (3, 1)))))

    _8_4_sub_ref = Set(map(x -> x[1], edges.(Ref(dgraph), ((8, 7), (8, 5), (5, 4), (7, 4)))))

    subs = extract_all!(compute_factorable_subgraphs(dgraph))
    _1_4_sub = subs[findfirst(x -> vertices(x) == (1, 4), subs)]
    _1_7_sub = subs[findfirst(x -> vertices(x) == (1, 7), subs)]
    _8_4_sub = subs[findfirst(x -> vertices(x) == (8, 4), subs)]
    _8_1_sub = subs[findfirst(x -> vertices(x) == (8, 1), subs)]
    _1_8_sub = subs[findfirst(x -> vertices(x) == (1, 8), subs)]

    @test issetequal(_1_4_sub_ref, subgraph_edges(_1_4_sub))
    factor_subgraph!(_1_4_sub)
    _1_7_sub_ref = Set(map(x -> x[1], edges.(Ref(dgraph), ((4, 1), (3, 1), (7, 4), (7, 6), (6, 3)))))


    @test issetequal(_1_7_sub_ref, subgraph_edges(_1_7_sub))
    @test issetequal(_8_4_sub_ref, subgraph_edges(_8_4_sub))
    factor_subgraph!(_8_4_sub)
    _8_1_sub_ref = Set(map(x -> x[1], edges.(Ref(dgraph), ((8, 7), (8, 4), (4, 1), (3, 1), (6, 3), (7, 6)))))
    @test issetequal(_8_1_sub_ref, subgraph_edges(_8_1_sub))
    @test issetequal(_8_1_sub_ref, subgraph_edges(_1_8_sub))

end

@testitem "subgraph_edges with branching" begin

    using FastDifferentiation.FDTests
    using FastDifferentiation.FDInternals

    @variables x

    x = Node(x)
    gr = DerivativeGraph((cos(x) * cos(x)) + x)
    # Vis.draw_dot(gr)
    # Vis.draw_dot(gr)
    sub = FactorableSubgraph{Int64,FD.DominatorSubgraph}(gr, 4, 1, BitVector([1]), BitVector([1]), BitVector([1]))

    edges_4_1 = collect(subgraph_edges(sub))

    sub = FactorableSubgraph{Int64,FD.PostDominatorSubgraph}(gr, 1, 4, BitVector([1]), BitVector([1]), BitVector([1]))
    edges_1_4 = collect(subgraph_edges(sub))

    @test count(x -> vertices(x) == (4, 3), edges_4_1) == 1
    @test count(x -> vertices(x) == (4, 1), edges_4_1) == 1
    @test count(x -> vertices(x) == (3, 2), edges_4_1) == 2
    @test count(x -> vertices(x) == (2, 1), edges_4_1) == 1

    @test count(x -> vertices(x) == (4, 3), edges_1_4) == 1
    @test count(x -> vertices(x) == (4, 1), edges_1_4) == 1
    @test count(x -> vertices(x) == (3, 2), edges_1_4) == 2
    @test count(x -> vertices(x) == (2, 1), edges_1_4) == 1
end

@testitem "deconstruct_subgraph" begin
    using FastDifferentiation.FDTests
    using FastDifferentiation.FDInternals

    graph, subs = simple_factorable_subgraphs()

    all_edges = collect(unique_edges(graph))

    _4_2 = all_edges[findfirst(x -> vertices(x) == (4, 2), all_edges)]
    _4_3 = all_edges[findfirst(x -> vertices(x) == (4, 3), all_edges)]
    _3_2 = all_edges[findfirst(x -> vertices(x) == (3, 2), all_edges)]
    _2_1 = all_edges[findfirst(x -> vertices(x) == (2, 1), all_edges)]
    _3_1 = all_edges[findfirst(x -> vertices(x) == (3, 1), all_edges)]

    ed, nod = deconstruct_subgraph(subs[1]) #can only deconstruct these two subgraphs because the larger ones need to be factored first.
    @test issetequal([4, 2, 3], nod)
    @test issetequal((_4_2, _4_3, _3_2), ed)

    ed, nod = deconstruct_subgraph(subs[2])
    @test issetequal((_3_2, _3_1, _2_1), ed)
    @test issetequal([1, 2, 3], nod)

    factor_subgraph!(subs[1]) #now can test larger subgraphs

    #new edges created and some edges deleted during factorization so get them again
    all_edges = collect(unique_edges(graph))

    _4_2 = all_edges[findfirst(x -> vertices(x) == (4, 2), all_edges)]
    _4_3 = all_edges[findfirst(x -> vertices(x) == (4, 3), all_edges)]
    _2_1 = all_edges[findfirst(x -> vertices(x) == (2, 1), all_edges)]
    _3_1 = all_edges[findfirst(x -> vertices(x) == (3, 1), all_edges)]

    ed, nod = deconstruct_subgraph(subs[3])
    println(ed)
    sub_4_1 = (_4_3, _4_2, _3_1, _2_1)
    @test issetequal(sub_4_1, ed)
    @test issetequal([1, 2, 3, 4], nod)
    ed, nod = deconstruct_subgraph(subs[4])
    @test issetequal(sub_4_1, ed)
    @test issetequal([1, 2, 3, 4], nod)
end

@testitem "subgraph reachable_roots, reachable_variables" begin

    using DataStructures
    using FastDifferentiation.FDInternals

    @variables nx1 ny2


    nxy3 = nx1 * ny2
    r2_4 = nx1 * nxy3
    r1_5 = r2_4 * nxy3

    gnodes = (nx1, ny2, nxy3, r2_4, r1_5)

    graph = DerivativeGraph([r1_5, r2_4])
    sub_heap = compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)

    subnums = ((5, 3), (4, 1), (5, 1), (1, 5), (3, 5), (1, 4))
    roots = (BitVector([1, 0]), BitVector([1, 1]), BitVector([1, 0]), BitVector([1, 0]), BitVector([1, 0]), BitVector([1, 1]))
    variables = (BitVector([1, 1]), BitVector([1, 0]), BitVector([1, 0]), BitVector([1, 0]), BitVector([1, 1]), BitVector([1, 0]))

    subgraphs = [x.subgraph for x in subs]
    #verify subgraphs have proper numbers in them
    for one_num in subnums
        @test one_num in subgraphs
    end

    for (i, one_root) in pairs(roots)
        sub = subs[findfirst(x -> x.subgraph == subnums[i], subs)]
        @test reachable_roots(sub) == one_root
    end

    for (i, one_var) in pairs(variables)
        sub = subs[findfirst(x -> x.subgraph == subnums[i], subs)]
        @test reachable_variables(sub) == one_var
    end
end

@testitem "Path_Iterator" begin

    using DataStructures
    using FastDifferentiation.FDInternals

    @variables nx1 ny2


    nxy3 = nx1 * ny2
    r2_4 = nx1 * nxy3
    r1_5 = r2_4 * nxy3

    gnodes = (nx1, ny2, nxy3, r2_4, r1_5)

    graph = DerivativeGraph([r1_5, r2_4])

    #first verify all nodes have the postorder numbers we expect
    for (i, nd) in pairs(gnodes)
        @test node(graph, i) == nd
    end

    sub_heap = compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)

    sub_5_3 = first(filter(x -> x.subgraph == (5, 3), subs))
    sub_3_5 = first((filter(x -> x.subgraph == (3, 5), subs)))

    rmask = reachable_dominance(sub_5_3)
    V = reachable_variables(sub_5_3)


    path_edges1 = [edges(graph, 4, 3)[1], edges(graph, 5, 4)[1]]
    path_edges2 = [edges(graph, 5, 3)[1]]


    start_edges = forward_edges(sub_5_3, dominated_node(sub_5_3))
    temp_edges = collect(edge_path(sub_5_3, start_edges[1]))

    @test all(x -> x[1] == x[2], zip(path_edges1, temp_edges))
    temp_edges = collect(edge_path(sub_5_3, start_edges[2]))
    @test all(x -> x[1] == x[2], zip(path_edges2, temp_edges))

    #for postdominator subgraph (3,5)

    start_edges = forward_edges(sub_3_5, dominated_node(sub_3_5))

    temp_edges = collect(edge_path(sub_3_5, start_edges[1]))
    @test all(x -> x[1] == x[2], zip(reverse(path_edges1), temp_edges))
    temp_edges = collect(edge_path(sub_3_5, start_edges[2]))
    @test all(x -> x[1] == x[2], zip(path_edges2, temp_edges))


    path_edges1 = [edges(graph, 4, 1)[1]]
    path_edges2 = [edges(graph, 3, 1)[1], edges(graph, 4, 3)[1]]
    sub_4_1 = first(filter(x -> x.subgraph == (4, 1), subs))
    sub_1_4 = first(filter(x -> x.subgraph == (1, 4), subs))


    start_edges = forward_edges(sub_4_1, dominated_node(sub_4_1))
    #for dominator subgraph (4,1)

    temp_edges = collect(edge_path(sub_4_1, start_edges[1]))
    @test all(x -> x[1] == x[2], zip(path_edges1, temp_edges))
    temp_edges = collect(edge_path(sub_4_1, start_edges[2]))
    @test all(x -> x[1] == x[2], zip(path_edges2, temp_edges))


    #for postdominator subgraph (1,4)
    start_edges = forward_edges(sub_1_4, dominated_node(sub_1_4))
    temp_edges = collect(edge_path(sub_1_4, start_edges[1]))

    @test all(x -> x[1] == x[2], zip(path_edges1, temp_edges))
    temp_edges = collect(edge_path(sub_1_4, start_edges[2]))
    @test all(x -> x[1] == x[2], zip(reverse(path_edges2), temp_edges))
end

@testitem "set_diff" begin
    using FastDifferentiation.FDInternals

    @test set_diff(falses(1), falses(1)) == falses(1)
    @test set_diff(falses(1), trues(1)) == falses(1)
    @test set_diff(trues(1), falses(1)) == trues(1)
    @test set_diff(trues(1), trues(1)) == falses(1)
end

@testitem "make_factored_edge" begin

    using DataStructures
    using FastDifferentiation.FDInternals

    @variables n1 n2


    n3 = n1 * n2
    n4 = n3 * n2
    n5 = n3 * n4
    n6 = n5 * n4

    graph = DerivativeGraph([n5, n6])

    sub_heap = compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)

    _5_3 = filter(x -> vertices(x) == (5, 3), subs)[1]
    e_5_3 = make_factored_edge(_5_3, evaluate_subgraph(_5_3))

    _3_5 = filter(x -> vertices(x) == (3, 5), subs)[1]
    e_3_5 = make_factored_edge(_3_5, evaluate_subgraph(_3_5))

    @test bit_equal(reachable_roots(e_5_3), BitVector([1, 0]))
    @test bit_equal(reachable_variables(e_5_3), BitVector([1, 1]))

    @test bit_equal(reachable_roots(e_3_5), BitVector([1, 1]))
    @test bit_equal(reachable_variables(e_3_5), BitVector([1, 0]))
end


@testitem "factor_subgraph simple ℝ²->ℝ²" begin

    using FastDifferentiation.FDInternals

    @variables x y

    nv1 = Node(x)
    nv2 = Node(y)
    n3 = nv1 * nv2
    n4 = n3 * nv1
    n5 = n3 * n4

    graph = DerivativeGraph([n4, n5])
    # factor_subgraph!(graph, postdominator_subgraph(2, 4, 2, BitVector([0, 1]), BitVector([0, 1])))
    subs = compute_factorable_subgraphs(graph)

    _5_3 = dominator_subgraph(graph, 5, 3, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _1_4 = postdominator_subgraph(graph, 1, 4, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _3_5 = postdominator_subgraph(graph, 3, 5, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _4_1 = dominator_subgraph(graph, 4, 1, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _5_1 = dominator_subgraph(graph, 5, 1, Bool[0, 1], Bool[0, 1], Bool[1, 0])
    _1_5 = postdominator_subgraph(graph, 1, 5, Bool[1, 0], Bool[0, 1], Bool[1, 0])

    sub_eval = evaluate_subgraph(_5_3)
    factor_subgraph!(_5_3)
end


@testitem "factor_subgraph 2" begin

    using FastDifferentiation.FDInternals

    @variables nx1 ny2


    n3 = nx1 * ny2
    n4 = n3 * ny2
    n5 = n3 * n4

    graph = DerivativeGraph([n5, n4])
    tmp = postdominator_subgraph(graph, 2, 4, BitVector([0, 1]), BitVector([0, 1]), BitVector([0, 1]))
    factor_subgraph!(tmp)
    @test length(edges(graph, 2, 4)) == 2

end



@testitem "evaluate_subgraph" begin
    using FastDifferentiation.FDTests
    using FastDifferentiation.FDInternals


    _, graph, _, _ = simple_dominator_graph()

    sub = postdominator_subgraph(graph, 1, 3, BitVector([1]), BitVector([1]), BitVector([1]))
end

@testitem "factor simple ℝ²->ℝ²" begin

    using FastDifferentiation.FDInternals

    @variables nx1 ny2


    n3 = nx1 * ny2
    r1_4 = sin(n3)
    r2_5 = cos(n3)

    graph = DerivativeGraph([r1_4, r2_5])
    result = _symbolic_jacobian!(graph, [nx1, ny2])

    #symbolic equality will work here because of common subexpression caching.
    @test result[1, 1] == cos(nx1 * ny2) * ny2
    @test result[1, 2] == cos(nx1 * ny2) * nx1
    @test result[2, 1] == -sin(nx1 * ny2) * ny2
    @test result[2, 2] == (-sin(nx1 * ny2)) * nx1
end


@testitem "subset" begin
    using FastDifferentiation.FDInternals

    a = falses(3)
    b = BitVector([1, 1, 1])
    @test subset(a, b)
    a = BitVector([1, 1, 1])
    @test subset(a, b)
    b = BitVector([1, 1, 0])
    @test !subset(a, b)
end


@testitem "constant and variable roots" begin

    using FastDifferentiation.FDInternals

    @variables x

    zr = Node(0.0)

    graph = DerivativeGraph([x, zr])
    jac = _symbolic_jacobian!(graph, [x])

    @test value(jac[1, 1]) == 1
    @test value(jac[2, 1]) == 0
end

@testitem "times_used PathEdge" begin
    using FastDifferentiation.FDInternals

    e = PathEdge(1, 2, Node(0), BitVector([1, 0, 1]), BitVector([0, 0, 1]))
    @test times_used(e) == 2
    e = PathEdge(1, 2, Node(0), BitVector([1, 0, 0]), BitVector([0, 0, 1]))
    @test times_used(e) == 1
    e = PathEdge(1, 2, Node(0), BitVector([1, 0, 1]), BitVector([1, 0, 1]))
    @test times_used(e) == 4
end

@testitem "path_sort_order" begin
    using FastDifferentiation.FDInternals
    e1 = PathEdge(1, 2, Node(0), BitVector([1, 0, 1]), BitVector([0, 0, 1]))
    e2 = PathEdge(3, 2, Node(0), BitVector([1, 0, 0]), BitVector([0, 0, 1]))
    @test path_sort_order(e1, e2) == true

    e3 = PathEdge(3, 2, Node(0), BitVector([1, 1, 0]), BitVector([0, 0, 1]))
    @test path_sort_order(e1, e3) == false
end

@testitem "multiply_sequence" begin

    using FastDifferentiation.FDInternals

    @variables x y z w u

    e1 = PathEdge(1, 2, x, BitVector([1, 0, 1]), BitVector([0, 0, 1]))
    e2 = PathEdge(3, 2, y, BitVector([1, 0, 0]), BitVector([0, 0, 1]))
    e3 = PathEdge(3, 2, z, BitVector([1, 1, 0]), BitVector([0, 0, 1]))
    e4 = PathEdge(3, 2, w, BitVector([1, 1, 0]), BitVector([1, 0, 1]))
    e5 = PathEdge(3, 2, u, BitVector([1, 1, 0]), BitVector([0, 1, 1]))


    path = [e1, e3, e2]   #2,2,1 times used
    @test (x * z) * y === multiply_sequence(path)
    path = [e4, e5]
    @test (w * u) === multiply_sequence(path)
    path = [e4, e5, e1]
    @test (w * u) * x === multiply_sequence(path)
    path = [e4, e5, e1, e3]
    @test (w * u) * (x * z) === multiply_sequence(path)
end


@testitem "factor ℝ¹->ℝ¹ " begin
    using FastDifferentiation.FDTests
    import FiniteDifferences
    using FastDifferentiation.FDInternals

    x, graph, _, _ = simple_dominator_graph()

    factor!(graph)
    fedge = edges(graph, 1, 4)[1]
    tmp0 = make_function([value(fedge)], [x])
    dfsimp(x) = tmp0([x])[1]
    x, graph, _, _ = simple_dominator_graph() #x is a new variable so have to make a new Node(x)

    tmp00 = make_function([root(graph, 1)], [x])
    origfsimp(x) = tmp00([x])[1]
    @test isapprox(FiniteDifferences.central_fdm(5, 1)(origfsimp, 3), dfsimp(3)[1])

    graph = complex_dominator_graph()
    factor!(graph)
    fedge = edges(graph, 1, 8)[1]
    tmp1 = make_function([value(fedge)], variables(graph))
    df(x) = tmp1(x)[1]

    graph = complex_dominator_graph()
    tmp2 = make_function([root(graph, 1)], variables(graph))
    origf(x) = tmp2(x)[1]

    for test_val in -3.0:0.013:3.0
        @test isapprox(FiniteDifferences.central_fdm(5, 1)(origf, test_val), df(test_val)[1])
    end
end


@testitem "jacobian" begin

    using FastDifferentiation.FDInternals

    @variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4

    graph = DerivativeGraph([n4, n5])

    df21(x, y) = 2 * x * y^3
    df22(x, y) = 4 * x^2 * y^2
    df11(x, y) = y^2
    df12(x, y) = 2 * x * y

    correct_jacobian = [df11 df12; df21 df22]
    copy_jac = _symbolic_jacobian(graph, [x, y])
    jac = _symbolic_jacobian!(graph, [x, y])

    @test all(copy_jac .== jac) #make sure the jacobian computed by copying the graph has the same variables as the one computed by destructively modifying the graph

    computed_jacobian = make_function(jac, [x, y])

    #verify the computed and hand caluclated jacobians agree.
    for _x in -1.0:0.01:1.0
        for _y in -1.0:0.3:1.0
            for index in CartesianIndices(correct_jacobian)
                @test isapprox(correct_jacobian[index](_x, _y), computed_jacobian([_x, _y])[index])
            end
        end
    end
end

@testitem "sparse_jacobian" begin
    using FastDifferentiation.FDTests
    import FastDifferentiation as FD

    @variables x y z

    sph_order = 10
    FD_graph = spherical_harmonics(sph_order, x, y, z)
    sprse = sparse_jacobian(FD.roots(FD_graph), [x, y, z])
    dense = jacobian(FD.roots(FD_graph), [x, y, z])

    for index in CartesianIndices(dense)
        if sprse[index] != dense[index] #empty elements in sprse get value Node{Int64,0} whereas zero elements in dense get value Node{Float64,0}. These are not == so need special case.
            @test FD.value(sprse[index]) == FD.value(dense[index])
        else
            @test sprse[index] == dense[index]
        end
    end
end

@testitem "sparse jacobian exe" begin
    using FastDifferentiation.FDTests
    using FastDifferentiation.FDInternals


    @variables x y z
    input_vars = [x, y, z]
    sph_order = 10
    FD_graph = spherical_harmonics(sph_order, x, y, z)
    sprse = sparse_jacobian(roots(FD_graph), input_vars)
    dense = jacobian(roots(FD_graph), input_vars)
    sprse_exe = make_function(sprse, input_vars)
    dense_exe = make_function(dense, input_vars)

    inputs = [1.0, 2.0, 3.0]

    sprse_res = sprse_exe(inputs)
    dense_res = dense_exe(inputs)

    @test isapprox(sprse_res, dense_res)
end

@testitem "spherical harmonics jacobian evaluation test" begin
    using FastDifferentiation.FDTests
    import FiniteDifferences
    using FastDifferentiation.FDInternals

    FD_graph = spherical_harmonics(10)
    mn_func = make_function(roots(FD_graph), variables(FD_graph))
    FD_func(vars...) = vec(mn_func(vars))

    graph_vars = variables(FD_graph)
    sym_func = make_function(FD.jacobian(roots(FD_graph), graph_vars), graph_vars)

    for xr in -1.0:0.3:1.0
        for yr in -1.0:0.3:1.0
            for zr = -1.0:0.3:1.0
                finite_diff = FiniteDifferences.jacobian(FiniteDifferences.central_fdm(12, 1, adapt=3), FD_func, xr, yr, zr)
                mat_form = hcat(finite_diff[1], finite_diff[2], finite_diff[3])
                symbolic = sym_func([xr, yr, zr])

                @test isapprox(symbolic, mat_form, rtol=1e-8)
            end
        end
    end
end

@testitem "Chebyshev jacobian evaluation test" begin
    import FiniteDifferences
    using FastDifferentiation.FDTests
    using FastDifferentiation.FDInternals

    chebyshev_order = 20
    FD_graph = chebyshev(FastSymbolic(), chebyshev_order)
    mn_func = make_function(roots(FD_graph), variables(FD_graph))
    FD_func(variables...) = vec(mn_func(variables...))

    func_wrap(x) = FD_func(x)[1]

    sym_func = make_function(FD.jacobian(roots(FD_graph), variables(FD_graph)), variables(FD_graph), in_place=false)

    for xr in -1.0:0.214:1.0
        finite_diff = FiniteDifferences.central_fdm(12, 1, adapt=3)(func_wrap, xr)

        symbolic = sym_func(xr)

        @test isapprox(symbolic[1, 1], finite_diff[1], rtol=1e-8)
    end

    tmp = Matrix{Float64}(undef, 1, 1)
    FD_graph = chebyshev(FastSymbolic(), chebyshev_order)
    sym_func = make_function(FD.jacobian(roots(FD_graph), variables(FD_graph)), variables(FD_graph), in_place=false)

    #the in place form of jacobian function
    for xr in -1.0:0.214:1.0
        finite_diff = FiniteDifferences.central_fdm(12, 1, adapt=3)(func_wrap, xr)

        symbolic = sym_func(xr, tmp)

        @test isapprox(symbolic[1, 1], finite_diff[1], rtol=1e-8)
    end
end

@testitem "derivative of matrix" begin

    using FastDifferentiation.FDInternals

    @variables nq1 nq2


    A = [
        cos(nq1) -cos(nq1)
        sin(nq1) sin(nq1)
    ]

    DA = [
        -sin(nq1) sin(nq1)
        cos(nq1) cos(nq1)
    ]

    @test isapprox(zeros(2, 2), value.(derivative(A, nq2))) #taking derivative wrt variable not present in the graph returns all zero matrix
    @test DA == derivative(A, nq1)
end

@testitem "jacobian_times_v" begin
    using FastDifferentiation.FDInternals
    using FastDifferentiation.FDTests

    order = 10

    FD_graph = spherical_harmonics(order)
    FD_func = roots(FD_graph)
    func_vars = variables(FD_graph)

    Jv, v_vars = jacobian_times_v(FD_func, func_vars)

    #compute the product the slow way
    Jv_slow = convert.(Node, jacobian(FD_func, func_vars) * v_vars)
    both_vars = [func_vars; v_vars]
    slow_symbolic = vec(reshape(Jv_slow, (length(Jv_slow), 1)))

    slow = make_function(slow_symbolic, both_vars)
    fast = make_function(Jv, both_vars)

    for _ in 1:100
        input = rand(length(func_vars) + length(v_vars))
        slow_val = slow(input)
        fast_val = fast(input)

        @test isapprox(slow_val, fast_val, rtol=1e-9)
    end

end

@testitem "jacobian_transpose_v" begin
    using FastDifferentiation.FDInternals
    using FastDifferentiation.FDTests

    order = 10

    FD_graph = spherical_harmonics(order)
    FD_func = roots(FD_graph)
    func_vars = variables(FD_graph)

    Jᵀv, r_vars = jacobian_transpose_v(FD_func, func_vars)

    Jᵀv_slow = convert.(Node, transpose(jacobian(FD_func, func_vars)) * r_vars)
    both_vars = [func_vars; r_vars]
    slow_symbolic = vec(reshape(Jᵀv_slow, (length(Jᵀv_slow), 1)))

    slow = make_function(slow_symbolic, both_vars)
    fast = make_function(Jᵀv, both_vars)

    for _ in 1:100
        input = rand(length(func_vars) + length(r_vars))
        slow_val = slow(input)
        fast_val = fast(input)

        @test isapprox(slow_val, fast_val, rtol=1e-8)
    end
end

@testitem "hessian" begin

    @variables x y z

    h = hessian(x^2 * y^2 * z^2, [x, y, z])
    h_exe = make_function(h, [x, y, z])

    @test isapprox(
        h_exe([1, 2, 3]),
        [
            72.0 72.0 48.0
            72.0 18.0 24.0
            48.0 24.0 8.0]
    )
end

@testitem "sparse hessian" begin
    using SparseArrays

    @variables x y z

    h = sparse_hessian(x^2 * y^2 * z^2, [x, y, z])
    h_exe = make_function(h, [x, y, z])
    inp = sparse(ones(Float64, 3, 3))

    @test isapprox(
        h_exe([1, 2, 3]),
        [
            72.0 72.0 48.0
            72.0 18.0 24.0
            48.0 24.0 8.0]
    )

    h2_exe = make_function(h, [x, y, z], in_place=true)
    h2_exe([1, 2, 3], inp)
    @test isapprox(inp,
        [
            72.0 72.0 48.0
            72.0 18.0 24.0
            48.0 24.0 8.0])
end

@testitem "hessian_times_v" begin
    using StaticArrays

    @variables x y
    v_vec = make_variables(:v, 2)
    f = x^2 * y^2
    hv_slow = convert.(FastDifferentiation.Node, hessian(f, [x, y]) * v_vec)

    hv_slow_exe = make_function(hv_slow, [[x, y]; v_vec])

    hv_fast, v_vec2 = hessian_times_v(f, [x, y])

    hv_fast_exe = make_function(hv_fast, [[x, y]; v_vec2])

    n_rand = 10
    for xval in rand(n_rand)
        for yval in rand(n_rand)
            for vval1 in rand(n_rand)
                for vval2 in rand(n_rand)
                    @test isapprox(hv_fast_exe(SVector{4,Float64}(xval, yval, vval1, vval2)), hv_slow_exe(SVector{4,Float64}(xval, yval, vval1, vval2)))
                end
            end
        end
    end
end

@testitem "sparse runtime generated functions" begin
    using SparseArrays

    @variables a11 a12 a13 a21 a22 a23 a31 a32 a33

    vars = vec([a11 a12 a13 a21 a22 a23 a31 a32 a33])
    spmat = sparse([a11 a12 a13; a21 a22 a23; a31 a32 a33])
    f1 = make_function(spmat, vars)
    inputs = [1 2 3 4 5 6 7 8 9]
    correct = [1 2 3; 4 5 6; 7 8 9]
    inp = similar(sprand(3, 3, 1.0))
    f2 = make_function(spmat, vars, in_place=true)

    @test f1(inputs) == correct
    f2(inputs, inp)
    @test inp == correct
end

@testitem "SArray return" begin
    using StaticArrays
    using FastDifferentiation.FDInternals
    using FastDifferentiation.FDTests

    @variables x y
    j = jacobian([x^2 * y^2, cos(x + y), log(x / y)], [x, y])
    j_exe = make_function(j, [x, y])
    @test typeof(j_exe([1.0, 2.0])) <: Array
    j_exe2 = make_function(SArray{Tuple{3,2}}(j), [x, y])
    @test typeof(j_exe2([1.0, 2.0])) <: StaticArray

    test_vec = [1.1, 2.3, 3.1]
    sph_func = spherical_harmonics(4)
    sph_jac = jacobian(roots(sph_func), variables(sph_func))
    mn_func1 = make_function(sph_jac, variables(sph_func)) #return type of executable should be Array
    m, n = size(sph_jac)
    mn_func2 = make_function(SMatrix{m,n}(sph_jac), variables(sph_func)) #return type of executable should be StaticArray
    @test typeof(mn_func1(test_vec)) <: Array
    @test typeof(mn_func2(test_vec)) <: StaticArray

    @test isapprox(mn_func1(test_vec), mn_func2(test_vec))
    @test isapprox(mn_func1(SVector{3}(test_vec)), mn_func2(test_vec))
    @test isapprox(mn_func1(SVector{3}(test_vec)), mn_func2(SVector{3}(test_vec)))
    @test isapprox(mn_func1(test_vec), mn_func2(SVector{3}(test_vec)))
end


end #module

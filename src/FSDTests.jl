module FSDTests
using DiffRules
using TermInterface
import SymbolicUtils
import Symbolics
using StaticArrays
using Memoize
using Rotations

using ..FastSymbolicDifferentiation

using TestItems

include("../FSDBenchmark/src/Chebyshev.jl")
include("../FSDBenchmark/src/SphericalHarmonics.jl")
include("../FSDBenchmark/src/LagrangianDynamics.jl")


""" Utility function for working with symbolic expressions as Symbolics.jl defines them."""
function number_of_operations(symbolic_expr)
    if SymbolicUtils.istree(symbolic_expr) && operation(symbolic_expr) ∈ (+, *, -)
        return 1 + sum(number_of_operations.(arguments(symbolic_expr)))
    else
        return 0
    end
end

function simple_dag(cache::Union{IdDict,Nothing}=IdDict())
    Symbolics.@variables zz y
    return expr_to_dag(zz^2 + y * (zz^2), cache), zz, y
end
export simple_dag

function simple_numbered_dag()
    Symbolics.@variables zz
    return expr_to_dag(zz * cos(zz^2))
end
export simple_numbered_dag

function dominators_dag()
    Symbolics.@variables zz
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
    Symbolics.@variables zz
    # generate node dag explicitly rather than starting from Symbolics expression because don't know how Symbolics might rearrange the nodes, causing the postorder numbers to change for tests.
    nx = FastSymbolicDifferentiation.Node(zz.val)
    sinx = FastSymbolicDifferentiation.Node(sin, MVector(nx))
    cosx = FastSymbolicDifferentiation.Node(cos, MVector(nx))
    A = FastSymbolicDifferentiation.Node(*, MVector(cosx, sinx))
    sinA = FastSymbolicDifferentiation.Node(sin, MVector(A))
    expsin = FastSymbolicDifferentiation.Node(*, MVector(A, FastSymbolicDifferentiation.Node(exp, MVector(sinx))))
    plus = FastSymbolicDifferentiation.Node(+, MVector(sinA, expsin))
    return plus
end
export complex_dominator_dag

complex_dominator_graph() = DerivativeGraph(complex_dominator_dag())
export complex_dominator_graph

function R2_R2_function()
    Symbolics.@variables x y

    nx = Node(x)
    ny = Node(y)
    n2 = nx * ny
    n4 = n2 * ny
    n5 = n2 * n4

    return DerivativeGraph([n5, n4])
end
export R2_R2_function

function simple_dominator_graph()
    Symbolics.@variables x

    nx = Node(x)
    ncos = Node(cos, nx)
    nplus = Node(+, ncos, nx)
    ntimes = Node(*, ncos, nplus)
    four_2_subgraph = Node(+, nplus, ncos)
    one_3_subgraph = Node(+, Node(*, Node(-1), Node(sin, nx)), Node(1))
    return x, DerivativeGraph(ntimes), four_2_subgraph, one_3_subgraph
end
export simple_dominator_graph

function simple_dominator_dgraph()
    x, graph, four_2_subgraph, one_3_subgraph = simple_dominator_graph()

    return graph, four_2_subgraph, one_3_subgraph, x
end
export simple_dominator_dgraph

@testitem "all_nodes" begin
    using FastSymbolicDifferentiation.FSDTests

    cache = IdDict()
    dag, x, y = simple_dag(cache)

    correct = expr_to_dag.([((x^2) + (y * (x^2))), (x^2), x, 2, (y * (x^2)), y], Ref(cache))
    tmp = all_nodes(dag)

    #verify that all the node expressions exist in the dag. Can't rely on them being in a particular order because Symbolics can
    #arbitrarily choose how to reorder trees.
    for expr in correct
        @test in(expr, tmp)
    end
end

# #TODO turn this into a real test
@testitem "make_function for Node" begin
    import Symbolics

    Symbolics.@variables x y

    A = [x^2+y 0 2x
        0 0 2y
        y^2+x 0 0]
    dag = expr_to_dag.(A)
    symbolics_answer = Symbolics.substitute.(A, Ref(Dict(x => 1.1, y => 2.3)))
    float_answer = similar(symbolics_answer, Float64)
    for index in eachindex(symbolics_answer)
        float_answer[index] = symbolics_answer[index].val
    end

    FSD_func = make_function.(dag, Ref([x, y]))
    res = [FSD_func[1, 1](1.1, 2.3) FSD_func[1, 2](0.0, 0.0) FSD_func[1, 3](1.1, 0.0)
        FSD_func[2, 1](0.0, 0.0) FSD_func[2, 2](0.0, 0.0) FSD_func[2, 3](0, 2.3)
        FSD_func[3, 1](1.1, 2.3) FSD_func[3, 2](0, 0) FSD_func[3, 3](0, 0)
    ]
    @test isapprox(res, float_answer)
end

@testitem "conversion from graph of FastSymbolicDifferentiation.Node to Symbolics expression" begin
    using FastSymbolicDifferentiation.FSDTests
    import Symbolics

    order = 7
    Symbolics.@variables x y z

    derivs = Symbolics.jacobian(SHFunctions(order, x, y, z), [x, y, z]; simplify=true)
    # show(@time SHDerivatives(order,x,y,z))
    tmp = expr_to_dag.(derivs)
    # show(@time expr_to_dag.(derivs))
    from_dag = dag_to_Symbolics_expression.(tmp)
    subs = Dict([x => rand(), y => rand(), z => rand()])
    @test isapprox(map(xx -> xx.val, Symbolics.substitute.(derivs, Ref(subs))), map(xx -> xx.val, Symbolics.substitute.(from_dag, Ref(subs))), atol=1e-14)
end

@testitem "is_tree" begin
    using Symbolics
    @variables x

    nx = Node(x)
    z = Node(0)
    tm = nx * nx

    @test is_tree(nx) == false
    @test is_leaf(nx) == true
    @test is_variable(nx) == true
    @test is_constant(nx) == false

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
    using Symbolics
    @variables x, y

    nx = Node(x)
    ny = Node(y)
    a = nx * ny
    @test derivative(a, Val(1)) == ny
    @test derivative(a, Val(2)) == nx
    @test derivative(nx) == Node(1)
    @test derivative(Node(1)) == Node(0)
end

@testitem "simple make_function" begin
    using FastSymbolicDifferentiation.FSDTests

    graph, four_2_subgraph, one_3_subgraph, _ = simple_dominator_dgraph()
    #not practical to compare the graphs directly since the order in which nodes come out of the differentiation
    #process is complicated. For dom subgraphs it depends on the order nodes appear in the parents list of a node. This 
    #is determined by code that has nothing to do with differentiation so don't want to take a dependency on it since it is
    #subject to change. Comparing graphs
    #directly would make the test fragile. Instead verify that the functions evaluate to the same thing.

    for testval in -π:0.08345:π
        correct_4_2_value = make_function(four_2_subgraph)
        computed_4_2_value = make_function(four_2_subgraph)
        correct_1_3_value = make_function(one_3_subgraph)
        computed_1_3_value = make_function(one_3_subgraph)
        @test isapprox(correct_4_2_value(testval), computed_4_2_value(testval), atol=1e-14)
        @test isapprox(correct_1_3_value(testval), computed_1_3_value(testval), atol=1e-14)
    end
end

@testitem "compute_factorable_subgraphs test order" begin
    using Symbolics

    @variables x y

    nv1 = Node(x)
    nv2 = Node(y)
    n3 = nv1 * nv2
    n4 = n3 * nv1

    n5 = n3 * n4

    graph = DerivativeGraph([n4, n5])
    # factor_subgraph!(graph, postdominator_subgraph(2, 4, 2, BitVector([0, 1]), BitVector([0, 1])))
    subs, _ = compute_factorable_subgraphs(graph)

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
    #test last two
    @test (value_equal(_5_1, subs[5]) && value_equal(_1_5, subs[6])) || (value_equal(_1_5, subs[5]) && value_equal(5_1, subs[6]))
end

@testitem "compute_factorable_subgraphs" begin
    using FastSymbolicDifferentiation.FSDTests

    dgraph = DerivativeGraph(complex_dominator_dag())

    subs, subs_dict = compute_factorable_subgraphs(dgraph)

    for sub in subs
        @test subs_dict[subgraph_vertices(sub)][1] == sub
    end

    for (sub, index) in values(subs_dict)
        @test sub in subs
    end

    equal_subgraphs(x, y) = dominating_node(x) == dominating_node(y) && dominated_node(x) == dominated_node(y) && times_used(x) == times_used(y) && dominance_mask(x) == dominance_mask(y)


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

@testitem "make_function" begin #test generation of derivative functions
    import Symbolics
    using Symbolics: substitute
    using FastSymbolicDifferentiation.FSDTests

    Symbolics.@variables zz
    dom_expr = zz * (cos(zz) + sin(zz))

    symbol_result = substitute(dom_expr, zz => 3.7)

    exe = make_function(dom_expr)

    @test exe(3.7) ≈ symbol_result

    Symbolics.@variables x, y
    sym_expr = cos(log(x) + sin(y)) * zz
    symbol_result = substitute(sym_expr, Dict([zz => 3.2, x => 2.5, y => 7.0]))

    exe = make_function(sym_expr)

    sym_expr2 = cos(2.0 * x) - sqrt(y)
    symbol_result = substitute(sym_expr2, Dict([x => 5.0, y => 3.2]))
    exe2 = make_function(sym_expr2, [x, y])

    symbol_val = symbol_result.val
    @test symbol_val ≈ exe2(5.0, 3.2)

    sym_expr3 = sin(x^3 + y^0.3)
    symbol_result3 = substitute(sym_expr3, Dict([x => 7.0, y => 0.4]))
    exe3 = make_function(sym_expr3, [x, y])
    symbol_val3 = symbol_result3.val
    @test symbol_val3 ≈ exe3(7.0, 0.4)

    #test to ensure that common terms are not reevaluated.
    sym_expr4 = sin(cos(x)) * cos(cos(x))
    symbol_result4 = substitute(sym_expr4, Dict([x => 7.0]))
    exe4 = make_function(sym_expr4, [x])
    symbol_val4 = symbol_result4.val
    @test symbol_val4 ≈ exe4(7.0)
end

@testitem "edges" begin
    using Symbolics

    @variables x y

    nx = Node(x)
    ny = Node(y)
    n2 = nx * ny
    n4 = n2 * ny
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
    using Symbolics

    @variables x, y

    #f = [cos(x) * sin(y), cos(x) + sin(y)]
    nx = Node(x)
    ny = Node(y)
    cosx = Node(cos, nx)
    sinx = Node(sin, nx)
    partial_cosx = Node(-, sinx)
    siny = Node(sin, ny)
    partial_siny = Node(cos, ny)
    ctimess = cosx * siny
    partial_times = [siny, cosx]
    cpluss = cosx + siny
    partial_plus = [Node(1), Node(1)]
    roots = [ctimess, cpluss]
    grnodes = [nx, ny, cosx, siny, cpluss, ctimess]

    correct_postorder_numbers = Dict((nx => 1, cosx => 2, ny => 3, siny => 4, ctimess => 5, cpluss => 6))

    graph = DerivativeGraph(roots)

    @test all([correct_postorder_numbers[node] == postorder_number(graph, node) for node in grnodes])

    correct_partials = Dict((cosx => [partial_cosx], siny => [partial_siny], ctimess => partial_times, cpluss => partial_plus))
    for (node, partials) in pairs(correct_partials)
        for (i, one_partial) in pairs(partials)
            f1 = dag_to_Symbolics_expression(partial_value(graph, node, i))
            f2 = dag_to_Symbolics_expression(one_partial)

            for test_point in BigFloat(-1):BigFloat(0.01):BigFloat(1) #graphs might have equivalent but different forms so evaluate at many points at high precision to verify derivatives are the same.
                v1 = Symbolics.value(substitute(f1, Dict((x => test_point), (y => test_point))))
                v2 = Symbolics.value(substitute(f2, Dict((x => test_point), (y => test_point))))
                @test isapprox(v1, v2, atol=1e-50)
            end
        end
    end
end

@testitem "DerivativeGraph pathmasks" begin
    using Symbolics
    @variables x, y

    nx = Node(x) #postorder # 2
    ny = Node(y) #postorder # 3
    xy = Node(*, nx, ny) #postorder # 4
    n5 = Node(5) #postorder # 1
    f1 = Node(*, n5, xy) #postorder 5
    n3 = Node(3) #postorder # 6
    f2 = Node(*, xy, n3) #postorder # 7
    roots = [f1, f2]
    graph = DerivativeGraph(roots)

    #nx,ny,xy shared by both roots. n5,f1 only on f1 path, n3,f2 only on f2 path
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
    using Symbolics
    @variables x, y

    #ℝ²->ℝ² function (f1,f2) = (5*(x*y),(x*y)*3)
    nx = Node(x) #postorder # 2
    ny = Node(y) #postorder # 3
    xy = Node(*, nx, ny) #postorder # 4
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
    for parent in relation_node_indices(piterator, postorder_number(graph, nx))
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
    @test children[1] == postorder_number(graph, nx)

    viterator = DomPathConstraint(graph, false, 2)
    children = Int64[]
    for child in relation_node_indices(viterator, postorder_number(graph, xy))
        push!(children, child)
    end
    @test length(children) == 1
    @test children[1] == 3
end

@testitem "edge_exists" begin
    using Symbolics
    @variables x, y

    nx = Node(x) #postorder # 2
    ny = Node(y) #postorder # 3
    xy = Node(*, nx, ny) #postorder # 4
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
    using Symbolics
    @variables x, y

    nx = Node(x) #postorder # 2
    ny = Node(y) #postorder # 3
    xy = Node(*, nx, ny) #postorder # 4
    n5 = Node(5) #postorder # 1
    f1 = Node(*, n5, xy) #postorder 5
    n3 = Node(3) #postorder # 6
    f2 = Node(*, xy, n3) #postorder # 7
    roots = [f1, f2]
    variables = [nx, ny]
    graph = DerivativeGraph(roots)

    previous_edges = unique_edges(graph)
    new_edge = PathEdge(1, 7, Node(y), length(variables), length(roots))
    FastSymbolicDifferentiation.add_edge!(graph, new_edge)

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

    #TODO need to add a test that deletes all the edges incident on a vertex and ensures that vertex is deleted.

    using Symbolics
    @variables x, y

    function reset_test(all_edges, graph, func::Function)
        for edge in all_edges
            tmp = func(edge)
            for i in eachindex(tmp)
                tmp[i] = 0
            end

            #can't delete edge till roots or variables are all unreachable

            FastSymbolicDifferentiation.delete_edge!(graph, edge)
            @test !edge_exists(graph, edge) #make sure edge has been deleted from graph

            delete!(all_edges, edge) #now delete edge and see if all the other edges that are still supposed to be in the graph are still there
            for edge2 in all_edges
                @test edge_exists(graph, edge2) #test other edges have not been deleted
            end
        end
    end
    nx = Node(x) #postorder # 2
    ny = Node(y) #postorder # 3
    xy = Node(*, nx, ny) #postorder # 4
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
    using Symbolics


    @variables x, y

    nx = Node(x) #postorder # 2
    ny = Node(y) #postorder # 3
    xy = Node(*, nx, ny) #postorder # 4
    n5 = Node(5) #postorder # 1
    f1 = Node(*, n5, xy) #postorder 5
    n3 = Node(3) #postorder # 6
    f2 = Node(*, xy, n3) #postorder # 7
    roots = [f1, f2]
    variables = [nx, ny]
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
    using Symbolics


    @variables x, y

    nx = Node(x) #postorder # 2
    ny = Node(y) #postorder # 3
    xy = Node(*, nx, ny) #postorder # 4
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


    nx = Node(x) #postorder #1
    ny = Node(y) #postorder #3
    ncos = Node(cos, nx) # 2
    nsin = Node(sin, ny) # 4
    n5 = Node(*, ncos, nsin) #5
    n6 = Node(*, n5, ny) #6
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
    using Symbolics

    using FastSymbolicDifferentiation: dom_subgraph, pdom_subgraph

    @variables x, y

    nx = Node(x) #postorder # 1
    ncos = Node(cos, nx) #pos
    ntimes1 = Node(*, ncos, nx)
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
    using Symbolics

    @variables x, y


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
    using FastSymbolicDifferentiation.FSDTests

    _, graph, four_2_subgraph, one_3_subgraph = simple_dominator_graph()



    subgraphs, subs_dict = compute_factorable_subgraphs(graph)
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

    @test index_1_3 < index_1_4
end

@testitem "TODO factor order with times used" begin
    #create test that checks order with different times used values for subgraphs
end


@testitem "subgraph reachable_roots, reachable_variables" begin
    using Symbolics

    @variables x, y

    nx1 = Node(x)
    ny2 = Node(y)
    nxy3 = nx1 * ny2
    r2_4 = nx1 * nxy3
    r1_5 = r2_4 * nxy3

    gnodes = (nx1, ny2, nxy3, r2_4, r1_5)

    graph = DerivativeGraph([r1_5, r2_4])
    subs, subs_dict = compute_factorable_subgraphs(graph)
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

@testitem "edges_on_path" begin
    using Symbolics

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
        @test node(graph, i) == nd
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
    @test all(x -> x[1] == x[2], zip(path_edges1, edges_on_path!(up_constraint, dominating, true, up_edges[1])))
    @test all(x -> x[1] == x[2], zip(path_edges2, edges_on_path!(up_constraint, dominating, true, up_edges[2])))

    #test for postdominator subgraph (3,5)
    down_constraint = PathConstraint(dominating_node(sub_3_5), graph, false, reachable_roots(sub_3_5), dominance_mask(sub_3_5))
    down_edges = relation_edges!(down_constraint, dominating)
    @test all(x -> x[1] == x[2], zip(reverse(path_edges1), edges_on_path!(down_constraint, dominated, false, down_edges[1])))
    @test all(x -> x[1] == x[2], zip(path_edges2, edges_on_path!(down_constraint, dominated, false, down_edges[2])))


    path_edges1 = [edges(graph, 4, 1)[1]]
    path_edges2 = [edges(graph, 3, 1)[1], edges(graph, 4, 3)[1]]
    sub_4_1 = first(filter(x -> x.subgraph == (4, 1), subs))
    sub_1_4 = first(filter(x -> x.subgraph == (1, 4), subs))
    dominating = sub_4_1.subgraph[1]
    dominated = sub_4_1.subgraph[2]

    #test for dominator subgraph (4,1)
    up_constraint = PathConstraint(dominating_node(sub_4_1), graph, true, dominance_mask(sub_4_1), reachable_variables(sub_4_1))
    up_edges = relation_edges!(up_constraint, dominated)
    @test all(x -> x[1] == x[2], zip(path_edges1, edges_on_path!(up_constraint, dominating, true, up_edges[1])))
    @test all(x -> x[1] == x[2], zip(path_edges2, edges_on_path!(up_constraint, dominating, true, up_edges[2])))

    dominating = sub_1_4.subgraph[1]
    dominated = sub_1_4.subgraph[2]
    #test for postdominator subgraph (1,4)
    down_constraint = PathConstraint(dominating_node(sub_1_4), graph, false, reachable_roots(sub_1_4), dominance_mask(sub_1_4))
    down_edges = relation_edges!(down_constraint, dominated)
    @test all(x -> x[1] == x[2], zip(path_edges1, edges_on_path!(down_constraint, dominating, false, down_edges[1])))
    @test all(x -> x[1] == x[2], zip(reverse(path_edges2), edges_on_path!(down_constraint, dominating, false, down_edges[2])))
end

@testitem "set_diff" begin
    @test set_diff(falses(1), falses(1)) == falses(1)
    @test set_diff(falses(1), trues(1)) == falses(1)
    @test set_diff(trues(1), falses(1)) == trues(1)
    @test set_diff(trues(1), trues(1)) == falses(1)
end

@testitem "edges of subgraph" begin
    using Symbolics

    function edges(a::FactorableSubgraph{T}, is_dominator::Bool) where {T}
        result = PathEdge{T}[]
        path_constraint = FastSymbolicDifferentiation.next_edge_constraint(a)

        for start_edge in relation_edges!(path_constraint, dominated_node(a))
            pedges = edges_on_path!(path_constraint, dominating_node(a), is_dominator, start_edge)
            append!(result, pedges)
        end
        return result
    end

    "returns edges in the subgraph `a` satisfying the constraint `dominance_mask(a) ⊆ reachable_roots(a)`"
    edges(a::FactorableSubgraph{T,FastSymbolicDifferentiation.DominatorSubgraph}) where {T} = edges(a, true)
    edges(a::FactorableSubgraph{T,FastSymbolicDifferentiation.PostDominatorSubgraph}) where {T} = edges(a, false)

    @variables x y

    nv1 = Node(x)
    nv2 = Node(y)
    n3 = nv1 * nv2
    n4 = n3 * nv1
    n5 = n3 * n4

    graph = DerivativeGraph([n4, n5])

    subs = compute_factorable_subgraphs(graph)[1]
    _5_3 = subs[findfirst(x -> x.subgraph == (5, 3), subs)]
    _4_1 = subs[findfirst(x -> x.subgraph == (4, 1), subs)]
    _1_4 = subs[findfirst(x -> x.subgraph == (1, 4), subs)]

    _5_3_edges = ((5, 3), (4, 3), (5, 4))
    _4_1_edges = ((4, 1), (4, 3), (3, 1))
    _1_4_edges = ((4, 1), (4, 3), (3, 1))

    @test length(edges(_5_3)) == 3
    @test length(edges(_4_1)) == 3
    @test length(edges(_1_4)) == 3

    function test_edges(subgraph, edge_list)
        for edge in edges(subgraph)
            @test (top_vertex(edge), bott_vertex(edge)) in edge_list
        end
    end

    for (subgraph, edges_of_subgraph) in zip((_5_3, _4_1, _1_4), (_5_3_edges, _4_1_edges), (_1_4_edges))
        test_edges(subgraph, edges_of_subgraph)
    end
end

@testitem "make_factored_edge" begin
    using Symbolics

    @variables v1, v2

    n1 = Node(v1)
    n2 = Node(v2)
    n3 = n1 * n2
    n4 = n3 * n2
    n5 = n3 * n4
    n6 = n5 * n4

    graph = DerivativeGraph([n5, n6])

    subs, sub_dict = compute_factorable_subgraphs(graph)

    _5_3 = filter(x -> subgraph_vertices(x) == (5, 3), subs)[1]
    e_5_3 = make_factored_edge(_5_3)

    _3_5 = filter(x -> subgraph_vertices(x) == (3, 5), subs)[1]
    e_3_5 = make_factored_edge(_3_5)

    @test bit_equal(reachable_roots(e_5_3), BitVector([1, 0]))
    @test bit_equal(reachable_variables(e_5_3), BitVector([1, 1]))

    @test bit_equal(reachable_roots(e_3_5), BitVector([1, 1]))
    @test bit_equal(reachable_variables(e_3_5), BitVector([1, 0]))
end


@testitem "factor_subgraph simple ℝ²->ℝ²" begin
    using Symbolics

    @variables x y

    nv1 = Node(x)
    nv2 = Node(y)
    n3 = nv1 * nv2
    n4 = n3 * nv1
    n5 = n3 * n4

    graph = DerivativeGraph([n4, n5])
    # factor_subgraph!(graph, postdominator_subgraph(2, 4, 2, BitVector([0, 1]), BitVector([0, 1])))
    subs, subgraph_dict = compute_factorable_subgraphs(graph)

    _5_3 = dominator_subgraph(graph, 5, 3, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _1_4 = postdominator_subgraph(graph, 1, 4, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _3_5 = postdominator_subgraph(graph, 3, 5, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _4_1 = dominator_subgraph(graph, 4, 1, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _5_1 = dominator_subgraph(graph, 5, 1, Bool[0, 1], Bool[0, 1], Bool[1, 0])
    _1_5 = postdominator_subgraph(graph, 1, 5, Bool[1, 0], Bool[0, 1], Bool[1, 0])

    sub_eval = evaluate_subgraph(_5_3)
    edges_to_delete, edge_to_add = factor_subgraph!(_5_3, sub_eval)
end


@testitem "factor_subgraph 2" begin
    using Symbolics
    @variables x y

    nx1 = Node(x)
    ny2 = Node(y)
    n3 = nx1 * ny2
    n4 = n3 * ny2
    n5 = n3 * n4

    graph = DerivativeGraph([n5, n4])
    tmp = postdominator_subgraph(graph, 2, 4, BitVector([0, 1]), BitVector([0, 1]), BitVector([0, 1]))
    factor_subgraph!(tmp, evaluate_subgraph(tmp))
    @test length(edges(graph, 2, 4)) == 1

end



@testitem "evaluate_subgraph" begin
    using FastSymbolicDifferentiation.FSDTests

    _, graph, _, _ = simple_dominator_graph()

    sub = postdominator_subgraph(graph, 1, 3, BitVector([1]), BitVector([1]), BitVector([1]))
end

@testitem "factor simple ℝ²->ℝ²" begin
    using Symbolics

    @variables x, y

    nx1 = Node(x)
    ny2 = Node(y)
    n3 = nx1 * ny2
    r1_4 = sin(n3)
    r2_5 = cos(n3)

    graph = DerivativeGraph([r1_4, r2_5])
    result = symbolic_jacobian!(graph, [nx1, ny2])

    @test result[1, 1] == cos(nx1 * ny2) * ny2
    @test result[1, 2] == cos(nx1 * ny2) * nx1
    @test result[2, 1] == -sin(nx1 * ny2) * ny2
    @test result[2, 2] == (-sin(nx1 * ny2)) * nx1
end

# @testitem "TODO: subgraph_exists" begin
#     @test false
# end

@testitem "subset" begin
    a = falses(3)
    b = BitVector([1, 1, 1])
    @test subset(a, b)
    a = BitVector([1, 1, 1])
    @test subset(a, b)
    b = BitVector([1, 1, 0])
    @test !subset(a, b)
end


@testitem "constants and variable roots" begin
    using Symbolics

    @variables x
    nx = Node(x)
    zr = Node(0.0)

    graph = DerivativeGraph([nx, zr])
    factor!(graph)
end

@testitem "times_used PathEdge" begin
    e = PathEdge(1, 2, Node(0), BitVector([1, 0, 1]), BitVector([0, 0, 1]))
    @test times_used(e) == 2
    e = PathEdge(1, 2, Node(0), BitVector([1, 0, 0]), BitVector([0, 0, 1]))
    @test times_used(e) == 1
    e = PathEdge(1, 2, Node(0), BitVector([1, 0, 1]), BitVector([1, 0, 1]))
    @test times_used(e) == 4
end

@testitem "TODO:subgraph_exists new version" begin end

@testitem "TODO:subgraph_nodes" begin end

@testitem "path_sort_order" begin
    e1 = PathEdge(1, 2, Node(0), BitVector([1, 0, 1]), BitVector([0, 0, 1]))
    e2 = PathEdge(3, 2, Node(0), BitVector([1, 0, 0]), BitVector([0, 0, 1]))
    @test path_sort_order(e1, e2) == true

    e3 = PathEdge(3, 2, Node(0), BitVector([1, 1, 0]), BitVector([0, 0, 1]))
    @test path_sort_order(e1, e3) == false
end

@testitem "multiply_sequence" begin
    using Symbolics
    @variables x, y, z, w, u

    e1 = PathEdge(1, 2, Node(x), BitVector([1, 0, 1]), BitVector([0, 0, 1]))
    e2 = PathEdge(3, 2, Node(y), BitVector([1, 0, 0]), BitVector([0, 0, 1]))
    e3 = PathEdge(3, 2, Node(z), BitVector([1, 1, 0]), BitVector([0, 0, 1]))
    e4 = PathEdge(3, 2, Node(w), BitVector([1, 1, 0]), BitVector([1, 0, 1]))
    e5 = PathEdge(3, 2, Node(u), BitVector([1, 1, 0]), BitVector([0, 1, 1]))


    path = [e1, e3, e2]   #2,2,1 times used
    @test (Node(x) * Node(z)) * Node(y) === multiply_sequence(path)
    path = [e4, e5]
    @test (Node(w) * Node(u)) === multiply_sequence(path)
    path = [e4, e5, e1]
    @test (Node(w) * Node(u)) * Node(x) === multiply_sequence(path)
    path = [e4, e5, e1, e3]
    @test (Node(w) * Node(u)) * (Node(x) * Node(z)) === multiply_sequence(path)
end


@testitem "factor ℝ¹->ℝ¹ " begin
    using FastSymbolicDifferentiation.FSDTests
    using FiniteDifferences

    _, graph, _, _ = simple_dominator_graph()
    factor!(graph)
    fedge = edges(graph, 1, 4)[1]
    dfsimp = make_function(edge_value(fedge))
    _, graph, _, _ = simple_dominator_graph()
    origfsimp = make_function(root(graph, 1))
    @test isapprox(central_fdm(5, 1)(origfsimp, 3), dfsimp(3))

    graph = complex_dominator_graph()
    factor!(graph)
    fedge = edges(graph, 1, 8)[1]
    df = make_function(edge_value(fedge))

    graph = complex_dominator_graph()
    origf = make_function(root(graph, 1))
    @test isapprox(central_fdm(5, 1)(origf, 3), df(3))
end

@testitem "factor ℝ²->ℝ² " begin
    using Symbolics

    @variables x y

    nx = Node(x)
    ny = Node(y)
    n2 = nx * ny
    n4 = n2 * ny
    n5 = n2 * n4

    graph = DerivativeGraph([n4, n5])
    # factor_subgraph!(graph, postdominator_subgraph(2, 4, 2, BitVector([0, 1]), BitVector([0, 1])))
    factor!(graph)
end

@testitem "symbolic_jacobian" begin
    using Symbolics

    @variables x y

    nx = Node(x)
    ny = Node(y)
    n2 = nx * ny
    n4 = n2 * ny
    n5 = n2 * n4

    graph = DerivativeGraph([n4, n5])

    df21(x, y) = 2 * x * y^3
    df22(x, y) = 3 * x^2 * y^2
    df11(x, y) = y^2
    df12(x, y) = 2 * x * y

    correct_jacobian = [df11 df12; df21 df22]

    computed_jacobian = make_function.(symbolic_jacobian!(graph, [nx, ny]), Ref([x, y]))


    #verify the computed and hand caluclated jacobians agree.
    for x in -1.0:0.01:1.0
        for y in -1.0:0.3:1.0
            for index in CartesianIndices(correct_jacobian)
                @test isapprox(correct_jacobian[index](x, y), computed_jacobian[index](x, y))
            end
        end
    end
end

@testitem "spherical harmonics jacobian evaluation test" begin
    using FastSymbolicDifferentiation.FSDTests
    using FiniteDifferences

    fsd_graph, x, y, z = to_graph(10)
    fsd_func = make_function(fsd_graph, Node.([x, y, z]))

    #hand computed derivative for order = 3
    # correct_derivatives(x, y, z) = [
    #     0.0 0.0 0.0
    #     0.0 0.0 0.0
    #     0.0 0.0 1.422074017360395
    #     0.0 0.0 0.0
    #     0.0 0.0 0.0
    #     (-0.3642161365141257*3*z) 0.0 (-0.3642161365141257*3*x)
    #     0.0 0.0 (2.0197963935867267*3*z)
    #     0.0 (-0.3642161365141257*3*z) (-0.3642161365141257*3*y)
    #     0.0 0.0 0.0
    # ]

    sym_func = jacobian_function!(fsd_graph, [Node(x), Node(y), Node(z)])


    for xr in -1.0:0.3:1.0
        for yr in -1.0:0.3:1.0
            for zr = -1.0:0.3:1.0
                finite_diff = jacobian(central_fdm(12, 1, adapt=3), fsd_func, xr, yr, zr)
                mat_form = hcat(finite_diff[1], finite_diff[2], finite_diff[3])
                symbolic = sym_func(xr, yr, zr)

                @test isapprox(symbolic, mat_form, rtol=1e-8)
            end
        end
    end
end

@testitem "Chebyshev jacobian evaluation test" begin
    using FiniteDifferences
    using FastSymbolicDifferentiation.FSDTests

    fsd_graph = chebyshev_graph(20)
    fsd_func = make_function(fsd_graph)

    func_wrap(x) = fsd_func(x)[1]

    sym_func = jacobian_function!(fsd_graph)

    for xr in -1.0:0.214:1.0
        finite_diff = central_fdm(12, 1, adapt=3)(func_wrap, xr)

        symbolic = sym_func(xr)

        @test isapprox(symbolic[1, 1], finite_diff[1], rtol=1e-8)
    end

    tmp = Matrix{Float64}(undef, 1, 1)
    fsd_graph = chebyshev_graph(20)
    sym_func = jacobian_function!(fsd_graph)

    #test the in place form of jacobian function
    for xr in -1.0:0.214:1.0
        finite_diff = central_fdm(12, 1, adapt=3)(func_wrap, xr)

        symbolic = sym_func(xr, tmp)

        @test isapprox(symbolic[1, 1], finite_diff[1], rtol=1e-8)
    end

end

@testitem "derivative of matrix" begin
    using Symbolics

    @variables q1 q2
    nq1 = Node(q1)
    nq2 = Node(q2)

    A = [
        cos(nq1) -cos(nq1)
        sin(nq1) sin(nq1)
    ]

    DA = [
        -sin(nq1) sin(nq1)
        cos(nq1) cos(nq1)
    ]

    @test isapprox(zeros(2, 2), node_value.(derivative(A, nq2))) #taking derivative wrt variable not present in the graph returns all zero matrix
    @test DA == derivative(A, nq1)
end

@testitem "derivative of Unspecified Function" begin
    using Symbolics
    using StaticArrays

    @variables x y

    ufn = function_of(:q, x, y)
    deriv = node_value(derivative(derivative(ufn, Val{1}()), Val{2}()))

    @assert deriv.variables == SVector(Node(x), Node(y))

    fn = DerivativeGraph(x * ufn)
    jac = symbolic_jacobian!(fn)
    @test jac[1, 1] == ufn + x * derivative(ufn, Val{1}())
    @test jac[1, 2] == x * derivative(ufn, Val{2}())
end

@testitem "unspecified function postorder" begin
    using Symbolics
    @variables x y

    q = function_of(:q, x, y)
    f = x * q + y * q
    graph = DerivativeGraph([f])

end


end #module
export Tests

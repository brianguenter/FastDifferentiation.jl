using TestItems

@testitem "FD.isa_connected_path 1" begin # case when path is one edge long
    using DataStructures
    import FastDifferentiation as FD


    FD.@variables x y

    func = x * x

    gr = FD.DerivativeGraph([func])
    subs_heap = FD.compute_factorable_subgraphs(gr)
    subs = extract_all!(subs_heap)
    test_sub = subs[1]

    etmp = FD.parent_edges(gr, FD.dominated_node(test_sub))
    rroots = FD.reachable_roots(etmp[1])
    rroots .= rroots .& .!rroots

    @test !FD.isa_connected_path(test_sub, etmp[1])
    @test FD.isa_connected_path(test_sub, etmp[2])
end

@testitem "FD.isa_connected_path 2" begin #cases when path is longer than one edge and various FD.edges have either FD.roots or FD.variables reset.

    using DataStructures
    import FastDifferentiation as FD

    FD.@variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4


    graph = FD.DerivativeGraph([n4, n5])
    subs_heap = FD.compute_factorable_subgraphs(graph)
    subs = extract_all!(subs_heap)
    println(subs)
    _5_3_index = findfirst(x -> FD.vertices(x) == (5, 3), subs)
    _5_3 = subs[_5_3_index]

    _2_4_index = findfirst(x -> FD.vertices(x) == (2, 4), subs)
    _2_4 = subs[_2_4_index]

    _3_5_index = findfirst(x -> FD.vertices(x) == (3, 5), subs)
    _3_5 = subs[_3_5_index]

    etmp = FD.edges(graph, 3, 5)[1]
    @test FD.isa_connected_path(_5_3, etmp)


    etmp = FD.edges(graph, 3, 4)[1]
    @test FD.isa_connected_path(_5_3, etmp)
    rts = FD.reachable_roots(etmp)
    rts[2] = 0

    @test !FD.isa_connected_path(_5_3, etmp)
    #reset path
    rts[2] = 1

    e2_4 = FD.edges(graph, 2, 4)[1]
    @test FD.isa_connected_path(_2_4, e2_4)
    e2_3 = FD.edges(graph, 2, 3)[1]
    @test FD.isa_connected_path(_2_4, e2_3)
    e3_4 = FD.edges(graph, 3, 4)[1]
    vars = FD.reachable_variables(e3_4)
    @. vars &= !vars
    @test !FD.isa_connected_path(_2_4, e3_4)
end
@testitem "add_non_dom_edges" begin
    import FastDifferentiation as FD
    using DataStructures


    #utility function to make it easier to create FD.edges and test them against FD.edges generated during graph operations.
    function edge_fields_equal(edge1, edge2)
        return edge1.top_vertex == edge2.top_vertex &&
               edge1.bott_vertex == edge2.bott_vertex &&
               edge1.edge_value == edge2.edge_value &&
               edge1.reachable_variables == edge2.reachable_variables &&
               edge1.reachable_roots == edge2.reachable_roots
    end

    FD.@variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4

    graph = FD.DerivativeGraph([n4, n5])
    subs_heap = FD.compute_factorable_subgraphs(graph)
    subs = extract_all!(subs_heap)
    _5_3 = subs[1]
    @test (5, 3) == FD.vertices(_5_3)

    FD.add_non_dom_edges!(_5_3)
    #single edge 3,4 should be split into two: ([r1,r2],[v1,v2]) -> ([r1],[v1,v2]),([r2],[v1,v2])
    edges3_4 = FD.edges(graph, 4, 3)
    @test length(edges3_4) == 2
    test_edge = FD.PathEdge(4, 3, y, BitVector([1, 1]), BitVector([0, 1]))
    @test count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1
    test_edge = (FD.PathEdge(4, 3, y, BitVector([1, 1]), BitVector([1, 0])))
    @test count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1

    graph = FD.DerivativeGraph([n4, n5])
    sub_heap = FD.compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)
    _2_4 = subs[2]
    @test (2, 4) == FD.vertices(_2_4)

    FD.add_non_dom_edges!(_2_4)
    #single edge 3,4 should be split in two: ([r1,r2],[v1,v2])->([r1,r2],[v1]),([r1,r2],[v2])
    edges3_4 = FD.edges(graph, 4, 3)
    @test length(edges3_4) == 2
    test_edge = FD.PathEdge(4, 3, y, BitVector([1, 0]), BitVector([1, 1]))
    @test count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1
    test_edge = (FD.PathEdge(4, 3, y, BitVector([0, 1]), BitVector([1, 1])))
    @test count(edge_fields_equal.(edges3_4, Ref(test_edge))) == 1
end

@testitem "iteration" begin
    import FastDifferentiation as FD
    using DataStructures


    FD.@variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4


    graph = FD.DerivativeGraph([n4, n5])
    subs_heap = FD.compute_factorable_subgraphs(graph)

    subs = extract_all!(subs_heap)

    _5_3_index = findfirst(x -> FD.vertices(x) == (5, 3), subs)
    _5_3 = subs[_5_3_index]

    _2_4_index = findfirst(x -> FD.vertices(x) == (2, 4), subs)
    _2_4 = subs[_2_4_index]

    _3_5_index = findfirst(x -> FD.vertices(x) == (3, 5), subs)
    _3_5 = subs[_3_5_index]

    e5_3 = FD.edges(graph, 5, 3)[1]

    pedges = collect(FD.edge_path(_5_3, e5_3))
    @test length(pedges) == 1
    @test e5_3 in pedges

    e3_4 = FD.edges(graph, 3, 4)[1]
    e5_4 = FD.edges(graph, 5, 4)[1]

    pedges = collect(FD.edge_path(_5_3, e3_4))
    @test length(pedges) == 2
    @test all(in.((e3_4, e5_4), Ref(pedges)))

    e2_3 = FD.edges(graph, 2, 3)[1]
    e2_4 = FD.edges(graph, 2, 4)[1]

    pedges = collect(FD.edge_path(_2_4, e3_4))
    @test length(pedges) == 2
    @test all(in.((e2_3, e3_4), Ref(pedges)))

    pedges = collect(FD.edge_path(_2_4, e2_4))
    @test length(pedges) == 1
    @test e2_4 in pedges
end



@testitem "FD.is_tree" begin
    import FastDifferentiation as FD

    FD.@variables x

    x = FD.Node(x)
    z = FD.Node(0)
    tm = x * x

    @test FD.is_tree(x) == false
    @test FD.is_leaf(x) == true
    @test FD.is_variable(x) == true
    @test FD.is_constant(x) == false

    @test FD.is_tree(z) == false
    @test FD.is_leaf(z) == true
    @test FD.is_constant(z) == true
    @test FD.is_variable(z) == false

    @test FD.is_leaf(tm) == false
    @test FD.is_variable(tm) == false
    @test FD.is_tree(tm) == true
    @test FD.is_constant(tm) == false
end

@testitem "derivative" begin
    import FastDifferentiation as FD
    FD.@variables x y


    a = x * y
    @test derivative(a, Val(1)) == y
    @test derivative(a, Val(2)) == x
end

@testitem "FD.compute_factorable_subgraphs test order" begin
    import FastDifferentiation as FD
    using DataStructures


    FD.@variables x y

    nv1 = FD.Node(x)
    nv2 = FD.Node(y)
    n3 = nv1 * nv2
    n4 = n3 * nv1

    n5 = n3 * n4

    graph = FD.DerivativeGraph([n4, n5])
    # FD.factor_subgraph!(graph, FD.postdominator_subgraph(2, 4, 2, BitVector([0, 1]), BitVector([0, 1])))
    sub_heap = FD.compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)

    _5_3 = FD.dominator_subgraph(graph, 5, 3, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _1_4 = FD.postdominator_subgraph(graph, 1, 4, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _3_5 = FD.postdominator_subgraph(graph, 3, 5, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _4_1 = FD.dominator_subgraph(graph, 4, 1, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _5_1 = FD.dominator_subgraph(graph, 5, 1, Bool[0, 1], Bool[0, 1], Bool[1, 0])
    _1_5 = FD.postdominator_subgraph(graph, 1, 5, Bool[1, 0], Bool[0, 1], Bool[1, 0])

    correctly_ordered_subs = (_5_3, _1_4, _3_5, _4_1, _5_1, _1_5) #order of last two could switch and still be correct but all others should be in exactly this order.

    tmp = zip(correctly_ordered_subs[1:4], subs[1:4])
    for (correct, computed) in tmp
        @test FD.value_equal(correct, computed)
    end
    #last two
    @test (FD.value_equal(_5_1, subs[5]) && FD.value_equal(_1_5, subs[6])) || (FD.value_equal(_1_5, subs[5]) && FD.value_equal(5_1, subs[6]))
end

@testitem "FD.compute_factorable_subgraphs" begin
    include("ShareTestCode.jl")
    using DataStructures

    import FastDifferentiation as FD
    dgraph = FD.DerivativeGraph(FDTests.complex_dominator_dag())

    sub_heap = FD.compute_factorable_subgraphs(dgraph)
    subs = extract_all!(sub_heap)


    equal_subgraphs(x, y) = FD.dominating_node(x) == FD.dominating_node(y) && FD.dominated_node(x) == FD.dominated_node(y) && FD.times_used(x) == FD.times_used(y) && FD.reachable_dominance(x) == FD.reachable_dominance(y)


    index_1_4 = findfirst(x -> equal_subgraphs(x, FD.postdominator_subgraph(dgraph, 1, 4, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    index_1_7 = findfirst(x -> equal_subgraphs(x, FD.postdominator_subgraph(dgraph, 1, 7, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    @test index_1_4 < index_1_7

    index_8_4 = findfirst(x -> equal_subgraphs(x, FD.dominator_subgraph(dgraph, 8, 4, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    index_8_3 = findfirst(x -> equal_subgraphs(x, FD.dominator_subgraph(dgraph, 8, 3, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    @test index_8_4 < index_8_3

    index_8_1d = findfirst(x -> equal_subgraphs(x, FD.dominator_subgraph(dgraph, 8, 1, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    @test index_8_4 < index_8_1d
    @test index_8_3 < index_8_1d
    @test index_1_7 < index_8_1d

    index_8_1p = findfirst(x -> equal_subgraphs(x, FD.dominator_subgraph(dgraph, 8, 1, BitVector([1]), BitVector([1]), BitVector([1]))), subs)
    @test index_8_4 < index_8_1p
    @test index_8_3 < index_8_1p
    @test index_1_7 < index_8_1p
end

@testitem "edges" begin
    import FastDifferentiation as FD


    FD.@variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4

    graph = FD.DerivativeGraph([n5, n4])

    function test_edge_access(graph, correct_num_edges, vert1, vert2)
        edge_group1 = FD.edges(graph, vert1, vert2)
        edge_group2 = FD.edges(graph, vert2, vert1)
        @test length(edge_group1) == correct_num_edges
        @test length(edge_group1) == length(edge_group2)
        for edge_pair in zip(edge_group1, edge_group2)
            @test edge_pair[1] == edge_pair[2]
        end

        for edge in edge_group1
            @test FD.top_vertex(edge) == max(vert1, vert2)
            @test FD.bott_vertex(edge) == min(vert1, vert2)
        end
    end

    edge_verts = ((1, 3), (3, 2), (2, 4), (5, 4), (3, 5), (4, 3))
    num_edges = (1, 1, 1, 1, 1, 1)

    for (edge, num) in zip(edge_verts, num_edges)
        test_edge_access(graph, num, edge[1], edge[2])
    end
end

@testitem "DerivativeGraph constructor" begin
    import FastDifferentiation as FD

    FD.@variables x y


    cosx = FD.Node(cos, x)
    sinx = FD.Node(sin, x)
    partial_cosx = FD.Node(-, sinx)
    siny = FD.Node(sin, y)
    partial_siny = FD.Node(cos, y)
    ctimess = cosx * siny
    partial_times = [siny, cosx]
    cpluss = cosx + siny
    partial_plus = [FD.Node(1), FD.Node(1)]
    rts = [ctimess, cpluss]
    grnodes = [x, y, cosx, siny, cpluss, ctimess]

    correct_postorder_numbers = Dict((x => 1, cosx => 2, y => 3, siny => 4, ctimess => 5, cpluss => 6))

    graph = FD.DerivativeGraph(rts)

    @test all([correct_postorder_numbers[node] == FD.postorder_number(graph, node) for node in grnodes])

    correct_partials = Dict((cosx => [partial_cosx], siny => [partial_siny], ctimess => partial_times, cpluss => partial_plus))
    #TODO need to use finite differences instead of Symbolics. 
    # for (node, partials) in pairs(correct_partials)
    #     for (i, one_partial) in pairs(partials)
    #         f1 = to_symbolics(partial_value(graph, node, i))
    #         f2 = to_symbolics(one_partial)

    #         for test_point in BigFloat(-1):BigFloat(0.01):BigFloat(1) #graphs might have equivalent but different forms so evaluate at many points at high precision to verify derivatives are the same.
    #             v1 = Symbolics.FD.value(Symbolics.substitute(f1, Dict((x => test_point), (y => test_point))))
    #             v2 = Symbolics.FD.value(Symbolics.substitute(f2, Dict((x => test_point), (y => test_point))))
    #             @test isapprox(v1, v2, atol=1e-50)
    #         end
    #     end
    # end
end

@testitem "DerivativeGraph pathmasks" begin

    import FastDifferentiation as FD

    FD.@variables x y

    xy = FD.Node(*, x, y) #postorder # 4
    n5 = FD.Node(5) #postorder # 1
    f1 = FD.Node(*, n5, xy) #postorder 5
    n3 = FD.Node(3) #postorder # 6
    f2 = FD.Node(*, xy, n3) #postorder # 7
    rts = [f1, f2]
    graph = FD.DerivativeGraph(rts)

    #x,y,xy shared by both FD.roots. n5,f1 only on f1 path, n3,f2 only on f2 path
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

    variable_path_masks = FD.compute_paths_to_variables(FD.num_vertices(graph), FD.edges(graph), FD.variable_index_to_postorder_number(graph))
    @test variable_path_masks == correct_variable_pathmasks
    parent_path_masks = FD.compute_paths_to_roots(FD.num_vertices(graph), FD.edges(graph), FD.root_index_to_postorder_number(graph))
    @test parent_path_masks == correct_roots_pathmasks
end

@testitem "ConstrainedPathIterator" begin
    import FastDifferentiation as FD

    FD.@variables x y

    #ℝ²->ℝ² function (f1,f2) = (5*(x*y),(x*y)*3)

    xy = FD.Node(*, x, y) #postorder # 4
    n5 = FD.Node(5) #postorder # 1
    f1 = FD.Node(*, n5, xy) #postorder 5
    n3 = FD.Node(3) #postorder # 6
    f2 = FD.Node(*, xy, n3) #postorder # 7
    rts = [f1, f2]
    graph = FD.DerivativeGraph(rts)
    root_masks = FD.compute_paths_to_roots(FD.num_vertices(graph), FD.edges(graph), FD.root_index_to_postorder_number(graph))
    variable_masks = FD.compute_paths_to_variables(FD.num_vertices(graph), FD.edges(graph), FD.variable_index_to_postorder_number(graph))

    piterator = FD.DomPathConstraint(graph, true, 1)

    parents = Int64[]
    correct_parents = (FD.postorder_number(graph, xy))
    for parent in FD.relation_node_indices(piterator, FD.postorder_number(graph, x))
        push!(parents, parent)
    end

    @test length(parents) == 1
    @test parents[1] == FD.postorder_number(graph, xy)

    piterator = FD.DomPathConstraint(graph, true, 1)
    parents = Int64[]
    pnum = FD.postorder_number(graph, xy)
    for parent in FD.relation_node_indices(piterator, pnum)
        push!(parents, parent)
    end
    @test length(parents) == 1
    @test FD.postorder_number(graph, f1) == parents[1]

    piterator = FD.DomPathConstraint(graph, true, 2)

    parents = Int64[]
    for parent in FD.relation_node_indices(piterator, FD.postorder_number(graph, xy))
        push!(parents, parent)
    end

    @test length(parents) == 1
    @test FD.postorder_number(graph, f2) == parents[1]


    viterator = FD.DomPathConstraint(graph, false, 1)
    children = Int64[]
    for child in FD.relation_node_indices(viterator, FD.postorder_number(graph, xy))
        push!(children, child)
    end
    @test length(children) == 1
    @test children[1] == FD.postorder_number(graph, x)

    viterator = FD.DomPathConstraint(graph, false, 2)
    children = Int64[]
    for child in FD.relation_node_indices(viterator, FD.postorder_number(graph, xy))
        push!(children, child)
    end
    @test length(children) == 1
    @test children[1] == 3
end

@testitem "edge_exists" begin
    import FastDifferentiation as FD


    FD.@variables x y

    xy = FD.Node(*, x, y) #postorder # 4
    n5 = FD.Node(5) #postorder # 1
    f1 = FD.Node(*, n5, xy) #postorder 5
    n3 = FD.Node(3) #postorder # 6
    f2 = FD.Node(*, xy, n3) #postorder # 7
    rts = [f1, f2]
    graph = FD.DerivativeGraph(rts)

    all_edges = FD.unique_edges(graph)

    for edge in all_edges
        @test FD.edge_exists(graph, edge)
    end

    @test !FD.edge_exists(graph, FD.PathEdge(1, 7, FD.Node(1), 2, 2)) #this edge is not in the graph

end

@testitem "add_edge! for FD.DerivativeGraph" begin
    import FastDifferentiation as FD


    FD.@variables x y

    xy = FD.Node(*, x, y) #postorder # 4
    n5 = FD.Node(5) #postorder # 1
    f1 = FD.Node(*, n5, xy) #postorder 5
    n3 = FD.Node(3) #postorder # 6
    f2 = FD.Node(*, xy, n3) #postorder # 7
    rts = [f1, f2]
    vars = [x, y]
    graph = FD.DerivativeGraph(rts)

    previous_edges = FD.unique_edges(graph)
    new_edge = FD.PathEdge(1, 7, FD.Node(y), length(vars), length(rts))
    FastDifferentiation.add_edge!(graph, new_edge)

    #make sure existing FD.edges are still in the graph.
    for edge in previous_edges
        @test FD.edge_exists(graph, edge)
    end

    prnts = FD.parents.(values(FD.edges(graph)))
    numprnts = sum(length.(prnts))
    chldrn = FD.children.(values(FD.edges(graph)))
    numchldrn = sum(length.(chldrn))
    num_edges = (numprnts + numchldrn) / 2
    @test num_edges == 7 #ensure number of FD.edges has increased by 1

    @test FD.edge_exists(graph, new_edge) #and that there is only one new edge
end

@testitem "delete_edge! for FD.DerivativeGraph" begin
    import FastDifferentiation as FD

    FD.@variables x y

    function reset_test(all_edges, graph, func::Function)
        for edge in all_edges
            tmp = func(edge)
            for i in eachindex(tmp)
                tmp[i] = 0
            end

            #can't delete edge till FD.roots or FD.variables are all unreachable

            FastDifferentiation.delete_edge!(graph, edge)
            @test !FD.edge_exists(graph, edge) #make sure edge has been deleted from graph

            delete!(all_edges, edge) #now delete edge and see if all the other FD.edges that are still supposed to be in the graph are still there
            for edge2 in all_edges
                @test FD.edge_exists(graph, edge2) #other FD.edges have not been deleted
            end
        end
    end

    xy = FD.Node(*, x, y) #postorder # 4
    n5 = FD.Node(5) #postorder # 1
    f1 = FD.Node(*, n5, xy) #postorder 5
    n3 = FD.Node(3) #postorder # 6
    f2 = FD.Node(*, xy, n3) #postorder # 7
    rts = [f1, f2]
    graph = FD.DerivativeGraph(rts)
    all_edges = FD.unique_edges(graph)

    reset_test(all_edges, graph, FD.reachable_roots)

    graph = FD.DerivativeGraph(rts)
    all_edges = FD.unique_edges(graph)

    reset_test(all_edges, graph, FD.reachable_variables)
end

@testitem "compute_edge_paths" begin

    import FastDifferentiation as FD


    FD.@variables x y

    xy = FD.Node(*, x, y) #postorder # 4
    n5 = FD.Node(5) #postorder # 1
    f1 = FD.Node(*, n5, xy) #postorder 5
    n3 = FD.Node(3) #postorder # 6
    f2 = FD.Node(*, xy, n3) #postorder # 7
    rts = [f1, f2]
    vars = [x, y]
    graph = FD.DerivativeGraph(rts)
    FD.compute_edge_paths!(FD.num_vertices(graph), FD.edges(graph), FD.variable_index_to_postorder_number(graph), FD.root_index_to_postorder_number(graph))

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

    for index in FD.each_vertex(graph)
        c_and_p = FD.node_edges(graph, index)
        for edge in [FD.parents(c_and_p); FD.children(c_and_p)]
            @test edge.reachable_variables == correct_variable_masks[(FD.top_vertex(edge), FD.bott_vertex(edge))]
            @test edge.reachable_roots == correct_root_masks[(FD.top_vertex(edge), FD.bott_vertex(edge))]
        end
    end
end

@testitem "dominators FD.DerivativeGraph" begin
    include("ShareTestCode.jl")
    import FastDifferentiation as FD



    FD.@variables x y

    xy = FD.Node(*, x, y) #postorder # 4
    n5 = FD.Node(5) #postorder # 1
    f1 = FD.Node(*, n5, xy) #postorder 5
    n3 = FD.Node(3) #postorder # 6
    f2 = FD.Node(*, xy, n3) #postorder # 7
    rts = [f1, f2]
    graph = FD.DerivativeGraph(rts)
    idoms = FDTests.compute_dominance_tables(graph, true)

    correct_dominators = [
        (1 => 5, 4 => 5, 2 => 4, 3 => 4, 5 => 5),
        (2 => 4, 3 => 4, 6 => 7, 4 => 7, 7 => 7)
    ]

    for (i, idom) in pairs(idoms)
        for elt in correct_dominators[i]
            @test elt[2] == idom[elt[1]]
        end
    end


    ncos = FD.Node(cos, x) # 2
    nsin = FD.Node(sin, y) # 4
    n5 = FD.Node(*, ncos, nsin) #5
    n6 = FD.Node(*, n5, y) #6
    n7 = FD.Node(*, n5, n6) #7
    nexp = FD.Node(exp, n6) # 8

    rts = [n7, nexp]
    graph = FD.DerivativeGraph(rts)
    idoms = FDTests.compute_dominance_tables(graph, true)

    correct_dominators = [
        (1 => 2, 2 => 5, 3 => 7, 4 => 5, 5 => 7, 6 => 7, 7 => 7),
        (1 => 2, 2 => 5, 3 => 6, 4 => 5, 5 => 6, 6 => 8, 8 => 8)
    ]

    for (i, idom) in pairs(idoms)
        for elt in correct_dominators[i]
            @test elt[2] == idom[elt[1]]
        end
    end

    pidoms = FDTests.compute_dominance_tables(graph, false)

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
    include("ShareTestCode.jl")
    import FastDifferentiation as FD

    using FastDifferentiation: dom_subgraph, pdom_subgraph


    FD.@variables x y

    ncos = FD.Node(cos, x) #pos
    ntimes1 = FD.Node(*, ncos, x)
    ntimes2 = FD.Node(*, ncos, ntimes1)
    rts = [ntimes2, ntimes1]
    graph = FD.DerivativeGraph(rts)
    idoms = FDTests.compute_dominance_tables(graph, true)
    pidoms = FDTests.compute_dominance_tables(graph, false)

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
    import FastDifferentiation as FD


    FD.@variables x y


    nx1 = FD.Node(x)
    ny2 = FD.Node(y)
    n3 = nx1 * ny2
    n4 = n3 * ny2
    n5 = n3 * n4

    graph = FD.DerivativeGraph([n5, n4])
    two = BitVector([1, 1])
    one_zero = BitVector([1, 0])
    zero_one = BitVector([0, 1])

    for node in (1, 2, 3, 4)
        @test FD.reachable_roots(graph, node) == two
    end
    @test FD.reachable_roots(graph, 5) == one_zero

    @test FD.reachable_variables(graph, 2) == zero_one
    for node in (3, 4, 5)
        @test FD.reachable_variables(graph, node) == two
    end

    @test FD.reachable_variables(graph, 1) == one_zero
end

@testitem "relation_edges" begin

end

@testitem "factor_order" begin
    include("ShareTestCode.jl")
    using DataStructures
    import FastDifferentiation as FD

    _, graph, four_2_subgraph, one_3_subgraph = FDTests.simple_dominator_graph()



    sub_heap = FD.compute_factorable_subgraphs(graph)
    subgraphs = extract_all!(sub_heap)
    @test length(subgraphs) == 4
    four_one = subgraphs[findfirst(x -> x.subgraph == (4, 1), subgraphs)]
    one_3 = subgraphs[findfirst(x -> x.subgraph == (1, 3), subgraphs)]
    four_2 = subgraphs[findfirst(x -> x.subgraph == (4, 2), subgraphs)]
    one_4 = subgraphs[findfirst(x -> x.subgraph == (1, 4), subgraphs)]

    @test FD.factor_order(one_3, four_one) == true
    @test FD.factor_order(four_2, four_one) == true
    @test FD.factor_order(one_3, four_2) == false
    @test FD.factor_order(four_2, one_3) == false
    @test FD.factor_order(four_one, one_3) == false
    @test FD.factor_order(four_one, four_2) == false


    @test FD.factor_order(one_3, one_4) == true
    @test FD.factor_order(four_2, one_4) == true
    @test FD.factor_order(one_3, four_2) == false
    @test FD.factor_order(four_2, one_3) == false
    @test FD.factor_order(one_4, one_3) == false
    @test FD.factor_order(one_4, four_2) == false #one_3 should be sorted before one_4

    equal_subgraphs(x, y) = FD.dominating_node(x) == FD.dominating_node(y) && FD.dominated_node(x) == FD.dominated_node(y) && FD.times_used(x) == FD.times_used(y) && FD.reachable_dominance(x) == FD.reachable_dominance(y)


    # doms = FD.dominator_subgraph.((
    #     (graph, 4, 2, BitVector([1, 0]), BitVector([1, 0])BitVector([1])),
    #     (graph, 4, 1, BitVector([1, 1]), BitVector([1, 0])BitVector([1]))))
    # pdoms = FD.postdominator_subgraph.((
    #     (graph, 1, 3, BitVector([1, 1]), BitVector([1]), BitVector([1, 1])),
    #     (graph, 1, 4, BitVector([1, 1]), BitVector([1]), BitVector([1, 0]))))
    # subs2 = collect((pdoms..., doms...))


    index_1_4 = findfirst(x -> equal_subgraphs(x, one_4), subgraphs)
    index_4_2 = findfirst(x -> equal_subgraphs(x, four_2), subgraphs)
    index_1_3 = findfirst(x -> equal_subgraphs(x, one_3), subgraphs)

    @test index_1_3 < index_1_4
end

@testitem "subgraph_edges" begin
    include("ShareTestCode.jl")
    using DataStructures
    import FastDifferentiation as FD

    dgraph = FD.DerivativeGraph([FDTests.complex_dominator_dag()])

    _1_4_sub_ref = Set(map(x -> x[1], FD.edges.(Ref(dgraph), ((4, 3), (4, 2), (2, 1), (3, 1)))))

    _8_4_sub_ref = Set(map(x -> x[1], FD.edges.(Ref(dgraph), ((8, 7), (8, 5), (5, 4), (7, 4)))))

    subs = extract_all!(FD.compute_factorable_subgraphs(dgraph))
    _1_4_sub = subs[findfirst(x -> FD.vertices(x) == (1, 4), subs)]
    _1_7_sub = subs[findfirst(x -> FD.vertices(x) == (1, 7), subs)]
    _8_4_sub = subs[findfirst(x -> FD.vertices(x) == (8, 4), subs)]
    _8_1_sub = subs[findfirst(x -> FD.vertices(x) == (8, 1), subs)]
    _1_8_sub = subs[findfirst(x -> FD.vertices(x) == (1, 8), subs)]

    @test issetequal(_1_4_sub_ref, FD.subgraph_edges(_1_4_sub))
    FD.factor_subgraph!(_1_4_sub)
    _1_7_sub_ref = Set(map(x -> x[1], FD.edges.(Ref(dgraph), ((4, 1), (3, 1), (7, 4), (7, 6), (6, 3)))))


    @test issetequal(_1_7_sub_ref, FD.subgraph_edges(_1_7_sub))
    @test issetequal(_8_4_sub_ref, FD.subgraph_edges(_8_4_sub))
    FD.factor_subgraph!(_8_4_sub)
    _8_1_sub_ref = Set(map(x -> x[1], FD.edges.(Ref(dgraph), ((8, 7), (8, 4), (4, 1), (3, 1), (6, 3), (7, 6)))))
    @test issetequal(_8_1_sub_ref, FD.subgraph_edges(_8_1_sub))
    @test issetequal(_8_1_sub_ref, FD.subgraph_edges(_1_8_sub))

end

@testitem "subgraph_edges with branching" begin
    import FastDifferentiation as FD
    include("ShareTestCode.jl")


    FD.@variables x

    x = FD.Node(x)
    gr = FD.DerivativeGraph((cos(x) * cos(x)) + x)
    # Vis.draw_dot(gr)
    # Vis.draw_dot(gr)
    sub = FD.FactorableSubgraph{Int64,FD.DominatorSubgraph}(gr, 4, 1, BitVector([1]), BitVector([1]), BitVector([1]))

    edges_4_1 = collect(FD.subgraph_edges(sub))

    sub = FD.FactorableSubgraph{Int64,FD.PostDominatorSubgraph}(gr, 1, 4, BitVector([1]), BitVector([1]), BitVector([1]))
    edges_1_4 = collect(FD.subgraph_edges(sub))

    @test count(x -> FD.vertices(x) == (4, 3), edges_4_1) == 1
    @test count(x -> FD.vertices(x) == (4, 1), edges_4_1) == 1
    @test count(x -> FD.vertices(x) == (3, 2), edges_4_1) == 2
    @test count(x -> FD.vertices(x) == (2, 1), edges_4_1) == 1

    @test count(x -> FD.vertices(x) == (4, 3), edges_1_4) == 1
    @test count(x -> FD.vertices(x) == (4, 1), edges_1_4) == 1
    @test count(x -> FD.vertices(x) == (3, 2), edges_1_4) == 2
    @test count(x -> FD.vertices(x) == (2, 1), edges_1_4) == 1
end

@testitem "deconstruct_subgraph" begin
    include("ShareTestCode.jl")
    import FastDifferentiation as FD

    graph, subs = FDTests.simple_factorable_subgraphs()

    all_edges = collect(FD.unique_edges(graph))

    _4_2 = all_edges[findfirst(x -> FD.vertices(x) == (4, 2), all_edges)]
    _4_3 = all_edges[findfirst(x -> FD.vertices(x) == (4, 3), all_edges)]
    _3_2 = all_edges[findfirst(x -> FD.vertices(x) == (3, 2), all_edges)]
    _2_1 = all_edges[findfirst(x -> FD.vertices(x) == (2, 1), all_edges)]
    _3_1 = all_edges[findfirst(x -> FD.vertices(x) == (3, 1), all_edges)]

    ed, nod = FD.deconstruct_subgraph(subs[1]) #can only deconstruct these two subgraphs because the larger ones need to be factored first.
    @test issetequal([4, 2, 3], nod)
    @test issetequal((_4_2, _4_3, _3_2), ed)

    ed, nod = FD.deconstruct_subgraph(subs[2])
    @test issetequal((_3_2, _3_1, _2_1), ed)
    @test issetequal([1, 2, 3], nod)

    FD.factor_subgraph!(subs[1]) #now can test larger subgraphs

    #new FD.edges created and some FD.edges deleted during factorization so get them again
    all_edges = collect(FD.unique_edges(graph))

    _4_2 = all_edges[findfirst(x -> FD.vertices(x) == (4, 2), all_edges)]
    _4_3 = all_edges[findfirst(x -> FD.vertices(x) == (4, 3), all_edges)]
    _2_1 = all_edges[findfirst(x -> FD.vertices(x) == (2, 1), all_edges)]
    _3_1 = all_edges[findfirst(x -> FD.vertices(x) == (3, 1), all_edges)]

    ed, nod = FD.deconstruct_subgraph(subs[3])
    println(ed)
    sub_4_1 = (_4_3, _4_2, _3_1, _2_1)
    @test issetequal(sub_4_1, ed)
    @test issetequal([1, 2, 3, 4], nod)
    ed, nod = FD.deconstruct_subgraph(subs[4])
    @test issetequal(sub_4_1, ed)
    @test issetequal([1, 2, 3, 4], nod)
end

@testitem "subgraph FD.reachable_roots, FD.reachable_variables" begin
    import FastDifferentiation as FD
    using DataStructures


    FD.@variables nx1 ny2


    nxy3 = nx1 * ny2
    r2_4 = nx1 * nxy3
    r1_5 = r2_4 * nxy3

    gnodes = (nx1, ny2, nxy3, r2_4, r1_5)

    graph = FD.DerivativeGraph([r1_5, r2_4])
    sub_heap = FD.compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)

    subnums = ((5, 3), (4, 1), (5, 1), (1, 5), (3, 5), (1, 4))
    rts = (BitVector([1, 0]), BitVector([1, 1]), BitVector([1, 0]), BitVector([1, 0]), BitVector([1, 0]), BitVector([1, 1]))
    vars = (BitVector([1, 1]), BitVector([1, 0]), BitVector([1, 0]), BitVector([1, 0]), BitVector([1, 1]), BitVector([1, 0]))

    subgraphs = [x.subgraph for x in subs]
    #verify subgraphs have proper numbers in them
    for one_num in subnums
        @test one_num in subgraphs
    end

    for (i, one_root) in pairs(rts)
        sub = subs[findfirst(x -> x.subgraph == subnums[i], subs)]
        @test FD.reachable_roots(sub) == one_root
    end

    for (i, one_var) in pairs(vars)
        sub = subs[findfirst(x -> x.subgraph == subnums[i], subs)]
        @test FD.reachable_variables(sub) == one_var
    end
end

@testitem "Path_Iterator" begin
    import FastDifferentiation as FD
    using DataStructures


    FD.@variables nx1 ny2


    nxy3 = nx1 * ny2
    r2_4 = nx1 * nxy3
    r1_5 = r2_4 * nxy3

    gnodes = (nx1, ny2, nxy3, r2_4, r1_5)

    graph = FD.DerivativeGraph([r1_5, r2_4])

    #first verify all nodes have the postorder numbers we expect
    for (i, nd) in pairs(gnodes)
        @test FD.node(graph, i) == nd
    end

    sub_heap = FD.compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)

    sub_5_3 = first(filter(x -> x.subgraph == (5, 3), subs))
    sub_3_5 = first((filter(x -> x.subgraph == (3, 5), subs)))

    rmask = FD.reachable_dominance(sub_5_3)
    V = FD.reachable_variables(sub_5_3)


    path_edges1 = [FD.edges(graph, 4, 3)[1], FD.edges(graph, 5, 4)[1]]
    path_edges2 = [FD.edges(graph, 5, 3)[1]]


    start_edges = FD.forward_edges(sub_5_3, FD.dominated_node(sub_5_3))
    temp_edges = collect(FD.edge_path(sub_5_3, start_edges[1]))

    @test all(x -> x[1] == x[2], zip(path_edges1, temp_edges))
    temp_edges = collect(FD.edge_path(sub_5_3, start_edges[2]))
    @test all(x -> x[1] == x[2], zip(path_edges2, temp_edges))

    #for postdominator subgraph (3,5)

    start_edges = FD.forward_edges(sub_3_5, FD.dominated_node(sub_3_5))

    temp_edges = collect(FD.edge_path(sub_3_5, start_edges[1]))
    @test all(x -> x[1] == x[2], zip(reverse(path_edges1), temp_edges))
    temp_edges = collect(FD.edge_path(sub_3_5, start_edges[2]))
    @test all(x -> x[1] == x[2], zip(path_edges2, temp_edges))


    path_edges1 = [FD.edges(graph, 4, 1)[1]]
    path_edges2 = [FD.edges(graph, 3, 1)[1], FD.edges(graph, 4, 3)[1]]
    sub_4_1 = first(filter(x -> x.subgraph == (4, 1), subs))
    sub_1_4 = first(filter(x -> x.subgraph == (1, 4), subs))


    start_edges = FD.forward_edges(sub_4_1, FD.dominated_node(sub_4_1))
    #for dominator subgraph (4,1)

    temp_edges = collect(FD.edge_path(sub_4_1, start_edges[1]))
    @test all(x -> x[1] == x[2], zip(path_edges1, temp_edges))
    temp_edges = collect(FD.edge_path(sub_4_1, start_edges[2]))
    @test all(x -> x[1] == x[2], zip(path_edges2, temp_edges))


    #for postdominator subgraph (1,4)
    start_edges = FD.forward_edges(sub_1_4, FD.dominated_node(sub_1_4))
    temp_edges = collect(FD.edge_path(sub_1_4, start_edges[1]))

    @test all(x -> x[1] == x[2], zip(path_edges1, temp_edges))
    temp_edges = collect(FD.edge_path(sub_1_4, start_edges[2]))
    @test all(x -> x[1] == x[2], zip(reverse(path_edges2), temp_edges))
end

@testitem "FD.set_diff" begin
    import FastDifferentiation as FD

    @test FD.set_diff(falses(1), falses(1)) == falses(1)
    @test FD.set_diff(falses(1), trues(1)) == falses(1)
    @test FD.set_diff(trues(1), falses(1)) == trues(1)
    @test FD.set_diff(trues(1), trues(1)) == falses(1)
end

@testitem "make_factored_edge" begin
    import FastDifferentiation as FD
    using DataStructures


    FD.@variables n1 n2


    n3 = n1 * n2
    n4 = n3 * n2
    n5 = n3 * n4
    n6 = n5 * n4

    graph = FD.DerivativeGraph([n5, n6])

    sub_heap = FD.compute_factorable_subgraphs(graph)
    subs = extract_all!(sub_heap)

    _5_3 = filter(x -> FD.vertices(x) == (5, 3), subs)[1]
    e_5_3 = FD.make_factored_edge(_5_3, FD.evaluate_subgraph(_5_3))

    _3_5 = filter(x -> FD.vertices(x) == (3, 5), subs)[1]
    e_3_5 = FD.make_factored_edge(_3_5, FD.evaluate_subgraph(_3_5))

    @test FD.bit_equal(FD.reachable_roots(e_5_3), BitVector([1, 0]))
    @test FD.bit_equal(FD.reachable_variables(e_5_3), BitVector([1, 1]))

    @test FD.bit_equal(FD.reachable_roots(e_3_5), BitVector([1, 1]))
    @test FD.bit_equal(FD.reachable_variables(e_3_5), BitVector([1, 0]))
end


@testitem "factor_subgraph simple ℝ²->ℝ²" begin
    import FastDifferentiation as FD


    FD.@variables x y

    nv1 = FD.Node(x)
    nv2 = FD.Node(y)
    n3 = nv1 * nv2
    n4 = n3 * nv1
    n5 = n3 * n4

    graph = FD.DerivativeGraph([n4, n5])
    # FD.factor_subgraph!(graph, FD.postdominator_subgraph(2, 4, 2, BitVector([0, 1]), BitVector([0, 1])))
    subs = FD.compute_factorable_subgraphs(graph)

    _5_3 = FD.dominator_subgraph(graph, 5, 3, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _1_4 = FD.postdominator_subgraph(graph, 1, 4, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _3_5 = FD.postdominator_subgraph(graph, 3, 5, Bool[0, 1], Bool[0, 1], Bool[1, 1])
    _4_1 = FD.dominator_subgraph(graph, 4, 1, Bool[1, 0], Bool[1, 1], Bool[1, 0])
    _5_1 = FD.dominator_subgraph(graph, 5, 1, Bool[0, 1], Bool[0, 1], Bool[1, 0])
    _1_5 = FD.postdominator_subgraph(graph, 1, 5, Bool[1, 0], Bool[0, 1], Bool[1, 0])

    sub_eval = FD.evaluate_subgraph(_5_3)
    FD.factor_subgraph!(_5_3)
end


@testitem "factor_subgraph 2" begin
    import FastDifferentiation as FD


    FD.@variables nx1 ny2


    n3 = nx1 * ny2
    n4 = n3 * ny2
    n5 = n3 * n4

    graph = FD.DerivativeGraph([n5, n4])
    tmp = FD.postdominator_subgraph(graph, 2, 4, BitVector([0, 1]), BitVector([0, 1]), BitVector([0, 1]))
    FD.factor_subgraph!(tmp)
    @test length(FD.edges(graph, 2, 4)) == 2

end



@testitem "evaluate_subgraph" begin
    include("ShareTestCode.jl")
    import FastDifferentiation as FD


    _, graph, _, _ = FDTests.simple_dominator_graph()

    sub = FD.postdominator_subgraph(graph, 1, 3, BitVector([1]), BitVector([1]), BitVector([1]))
end

@testitem "factor simple ℝ²->ℝ²" begin
    import FastDifferentiation as FD


    FD.@variables nx1 ny2


    n3 = nx1 * ny2
    r1_4 = sin(n3)
    r2_5 = cos(n3)

    graph = FD.DerivativeGraph([r1_4, r2_5])
    result = FD._symbolic_jacobian!(graph, [nx1, ny2])

    #symbolic equality will work here because of common subexpression caching.
    @test result[1, 1] == cos(nx1 * ny2) * ny2
    @test result[1, 2] == cos(nx1 * ny2) * nx1
    @test result[2, 1] == -sin(nx1 * ny2) * ny2
    @test result[2, 2] == (-sin(nx1 * ny2)) * nx1
end


@testitem "FD.subset" begin
    import FastDifferentiation as FD

    a = falses(3)
    b = BitVector([1, 1, 1])
    @test FD.subset(a, b)
    a = BitVector([1, 1, 1])
    @test FD.subset(a, b)
    b = BitVector([1, 1, 0])
    @test !FD.subset(a, b)
end


@testitem "constant and variable FD.roots" begin
    import FastDifferentiation as FD


    FD.@variables x

    zr = FD.Node(0.0)

    graph = FD.DerivativeGraph([x, zr])
    jac = FD._symbolic_jacobian!(graph, [x])

    @test FD.value(jac[1, 1]) == 1
    @test FD.value(jac[2, 1]) == 0
end

@testitem "FD.times_used FD.PathEdge" begin
    import FastDifferentiation as FD

    e = FD.PathEdge(1, 2, FD.Node(1), BitVector([1, 0, 1]), BitVector([0, 0, 1]))
    @test FD.times_used(e) == 2
    e = FD.PathEdge(1, 2, FD.Node(1), BitVector([1, 0, 0]), BitVector([0, 0, 1]))
    @test FD.times_used(e) == 1
    e = FD.PathEdge(1, 2, FD.Node(1), BitVector([1, 0, 1]), BitVector([1, 0, 1]))
    @test FD.times_used(e) == 4
end

@testitem "FD.path_sort_order" begin
    import FastDifferentiation as FD

    e1 = FD.PathEdge(1, 2, FD.Node(1), BitVector([1, 0, 1]), BitVector([0, 0, 1]))
    e2 = FD.PathEdge(3, 2, FD.Node(1), BitVector([1, 0, 0]), BitVector([0, 0, 1]))
    @test FD.path_sort_order(e1, e2) == true

    e3 = FD.PathEdge(3, 2, FD.Node(1), BitVector([1, 1, 0]), BitVector([0, 0, 1]))
    @test FD.path_sort_order(e1, e3) == false
end

@testitem "FD.multiply_sequence" begin
    import FastDifferentiation as FD


    FD.@variables x y z w u

    e1 = FD.PathEdge(1, 2, x, BitVector([1, 0, 1]), BitVector([0, 0, 1]))
    e2 = FD.PathEdge(3, 2, y, BitVector([1, 0, 0]), BitVector([0, 0, 1]))
    e3 = FD.PathEdge(3, 2, z, BitVector([1, 1, 0]), BitVector([0, 0, 1]))
    e4 = FD.PathEdge(3, 2, w, BitVector([1, 1, 0]), BitVector([1, 0, 1]))
    e5 = FD.PathEdge(3, 2, u, BitVector([1, 1, 0]), BitVector([0, 1, 1]))


    path = [e1, e3, e2]   #2,2,1 times used
    @test (x * z) * y === FD.multiply_sequence(path)
    path = [e4, e5]
    @test (w * u) === FD.multiply_sequence(path)
    path = [e4, e5, e1]
    @test (w * u) * x === FD.multiply_sequence(path)
    path = [e4, e5, e1, e3]
    @test (w * u) * (x * z) === FD.multiply_sequence(path)
end


@testitem "factor ℝ¹->ℝ¹ " begin
    include("ShareTestCode.jl")
    import FiniteDifferences
    import FastDifferentiation as FD

    x, graph, _, _ = FDTests.simple_dominator_graph()

    FD.factor!(graph)
    fedge = FD.edges(graph, 1, 4)[1]
    tmp0 = FD.make_function([FD.value(fedge)], [x])
    dfsimp(x) = tmp0([x])[1]
    x, graph, _, _ = FDTests.simple_dominator_graph() #x is a new variable so have to make a new FD.Node(x)

    tmp00 = FD.make_function([FD.root(graph, 1)], [x])
    origfsimp(x) = tmp00([x])[1]
    @test isapprox(FiniteDifferences.central_fdm(5, 1)(origfsimp, 3), dfsimp(3)[1])

    graph = FDTests.complex_dominator_graph()
    FD.factor!(graph)
    fedge = FD.edges(graph, 1, 8)[1]
    tmp1 = FD.make_function([FD.value(fedge)], FD.variables(graph))
    df(x) = tmp1(x)[1]

    graph = FDTests.complex_dominator_graph()
    tmp2 = FD.make_function([FD.root(graph, 1)], FD.variables(graph))
    origf(x) = tmp2(x)[1]

    for test_val in -3.0:0.013:3.0
        @test isapprox(FiniteDifferences.central_fdm(5, 1)(origf, test_val), df(test_val)[1])
    end
end


@testitem "jacobian" begin
    import FastDifferentiation as FD


    FD.@variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4

    graph = FD.DerivativeGraph([n4, n5])

    df21(x, y) = 2 * x * y^3
    df22(x, y) = 4 * x^2 * y^2
    df11(x, y) = y^2
    df12(x, y) = 2 * x * y

    correct_jacobian = [df11 df12; df21 df22]
    copy_jac = FD._symbolic_jacobian(graph, [x, y])
    jac = FD._symbolic_jacobian!(graph, [x, y])

    @test all(copy_jac .== jac) #make sure the jacobian computed by copying the graph has the same FD.variables as the one computed by destructively modifying the graph

    computed_jacobian = FD.make_function(jac, [x, y])

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
    include("ShareTestCode.jl")
    import FastDifferentiation as FD

    FD.@variables x y z

    sph_order = 10
    FD_graph = FDTests.spherical_harmonics(sph_order, x, y, z)
    sprse = sparse_jacobian(FD.roots(FD_graph), [x, y, z])
    dense = jacobian(FD.roots(FD_graph), [x, y, z])

    for index in CartesianIndices(dense)
        @test !xor(FD.is_zero(sprse[index]), FD.is_zero(dense[index])) #See https://vscode.dev/github/brianguenter/FastDifferentiation.jl/blob/main/src/ExpressionGraph.jl#L271. sprse zero elements
        # are FD.Node{Int64,0} while dense elements are FD.Node{Float64,0}.
        if !FD.is_zero(sprse[index])
            @test sprse[index] === dense[index]
        end
    end
end

@testitem "sparse jacobian exe" begin
    include("ShareTestCode.jl")
    import FastDifferentiation as FD


    FD.@variables x y z
    input_vars = [x, y, z]
    sph_order = 10
    FD_graph = FDTests.spherical_harmonics(sph_order, x, y, z)
    sprse = sparse_jacobian(FD.roots(FD_graph), input_vars)
    dense = jacobian(FD.roots(FD_graph), input_vars)
    sprse_exe = FD.make_function(sprse, input_vars)
    dense_exe = FD.make_function(dense, input_vars)

    inputs = [1.0, 2.0, 3.0]

    sprse_res = sprse_exe(inputs)
    dense_res = dense_exe(inputs)

    @test isapprox(sprse_res, dense_res)
end

@testitem "spherical harmonics jacobian evaluation test" begin
    include("ShareTestCode.jl")
    import FiniteDifferences
    import FastDifferentiation as FD

    FD_graph = FDTests.spherical_harmonics(10)
    mn_func = FD.make_function(FD.roots(FD_graph), FD.variables(FD_graph))
    FD_func(vars...) = vec(mn_func(vars))

    graph_vars = FD.variables(FD_graph)
    sym_func = FD.make_function(FD.jacobian(FD.roots(FD_graph), graph_vars), graph_vars)

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
    include("ShareTestCode.jl")
    import FastDifferentiation as FD

    chebyshev_order = 20
    FD_graph = FDTests.chebyshev(FDTests.FastSymbolic(), chebyshev_order)
    mn_func = FD.make_function(FD.roots(FD_graph), FD.variables(FD_graph))
    FD_func(variables...) = vec(mn_func(variables...))

    func_wrap(x) = FD_func(x)[1]

    sym_func = FD.make_function(FD.jacobian(FD.roots(FD_graph), FD.variables(FD_graph)), FD.variables(FD_graph), in_place=false)

    for xr in -1.0:0.214:1.0
        finite_diff = FiniteDifferences.central_fdm(12, 1, adapt=3)(func_wrap, xr)

        symbolic = sym_func(xr)

        @test isapprox(symbolic[1, 1], finite_diff[1], rtol=1e-8)
    end

    tmp = Matrix{Float64}(undef, 1, 1)
    FD_graph = FDTests.chebyshev(FDTests.FastSymbolic(), chebyshev_order)
    sym_func = FD.make_function(FD.jacobian(FD.roots(FD_graph), FD.variables(FD_graph)), FD.variables(FD_graph), in_place=false)

    #the in place form of jacobian function
    for xr in -1.0:0.214:1.0
        finite_diff = FiniteDifferences.central_fdm(12, 1, adapt=3)(func_wrap, xr)

        symbolic = sym_func(xr, tmp)

        @test isapprox(symbolic[1, 1], finite_diff[1], rtol=1e-8)
    end
end

@testitem "derivative of matrix" begin
    import FastDifferentiation as FD


    FD.@variables nq1 nq2


    A = [
        cos(nq1) -cos(nq1)
        sin(nq1) sin(nq1)
    ]

    DA = [
        -sin(nq1) sin(nq1)
        cos(nq1) cos(nq1)
    ]

    @test isapprox(zeros(2, 2), FD.value.(derivative(A, nq2))) #taking derivative wrt variable not present in the graph returns all zero matrix
    @test DA == derivative(A, nq1)
end

@testitem "jacobian_times_v" begin
    import FastDifferentiation as FD
    using Random
    include("ShareTestCode.jl")

    order = 10

    FD_graph = FDTests.spherical_harmonics(order)
    FD_func = FD.roots(FD_graph)
    func_vars = FD.variables(FD_graph)

    Jv, v_vars = FD.jacobian_times_v(FD_func, func_vars)

    #compute the product the slow way
    Jv_slow = convert.(FD.Node, jacobian(FD_func, func_vars) * v_vars)
    both_vars = [func_vars; v_vars]
    slow_symbolic = vec(reshape(Jv_slow, (length(Jv_slow), 1)))

    slow = FD.make_function(slow_symbolic, both_vars)
    fast = FD.make_function(Jv, both_vars)

    for _ in 1:100
        input = rand(length(func_vars) + length(v_vars))
        slow_val = slow(input)
        fast_val = fast(input)

        @test isapprox(slow_val, fast_val, rtol=1e-9)
    end

    #test for https://github.com/brianguenter/FastDifferentiation.jl/issues/60
    function heat(u, p, t)
        gamma, dx = p
        dxxu = (circshift(u, 1) .- 2 * u .+ circshift(u, -1)) / dx^2
        return gamma * dxxu
    end


    Ns = 5
    xmin, xmax = 0.0, 1.0
    dx = (xmax - xmin) / (Ns - 1)

    xs = collect(LinRange(xmin, xmax, Ns))
    using FastDifferentiation
    FastDifferentiation.@variables uv pv tv
    uvs = make_variables(:uv, Ns)
    pvs = make_variables(:pv, 2)
    tvs = make_variables(:tv, 1)
    fv = heat(uvs, pvs, tvs)
    js = FastDifferentiation.jacobian(fv, uvs) # The Jacobian js is of shape Ns x Ns

    vjps, rvs = FastDifferentiation.jacobian_times_v(fv, uvs) # The variable vjps is of shape Ns+2 and rvs is of shape Ns

    input_vec = [uvs; pvs; tvs; rvs]
    jv = make_function(FastDifferentiation.Node.(js * rvs), input_vec)
    jtv2 = make_function(vjps, input_vec)

    for i in 1:100
        rng = Random.Xoshiro(123)
        tvec = rand(rng, length(input_vec))
        @test isapprox(jv(tvec), jtv2(tvec))
    end
end

@testitem "jacobian_transpose_v" begin
    import FastDifferentiation as FD
    using Random
    include("ShareTestCode.jl")

    order = 10

    FD_graph = FDTests.spherical_harmonics(order)
    FD_func = FD.roots(FD_graph)
    func_vars = FD.variables(FD_graph)

    Jᵀv, r_vars = jacobian_transpose_v(FD_func, func_vars)

    Jᵀv_slow = convert.(FD.Node, transpose(jacobian(FD_func, func_vars)) * r_vars)
    both_vars = [func_vars; r_vars]
    slow_symbolic = vec(reshape(Jᵀv_slow, (length(Jᵀv_slow), 1)))

    slow = FD.make_function(slow_symbolic, both_vars)
    fast = FD.make_function(Jᵀv, both_vars)

    for _ in 1:100
        input = rand(length(func_vars) + length(r_vars))
        slow_val = slow(input)
        fast_val = fast(input)

        @test isapprox(slow_val, fast_val, rtol=1e-8)
    end

    #test for https://github.com/brianguenter/FastDifferentiation.jl/issues/60
    function heat(u, p, t)
        gamma, dx = p
        dxxu = (circshift(u, 1) .- 2 * u .+ circshift(u, -1)) / dx^2
        return gamma * dxxu
    end


    Ns = 5
    xmin, xmax = 0.0, 1.0
    dx = (xmax - xmin) / (Ns - 1)

    xs = collect(LinRange(xmin, xmax, Ns))
    using FastDifferentiation
    FastDifferentiation.@variables uv pv tv
    uvs = make_variables(:uv, Ns)
    pvs = make_variables(:pv, 2)
    tvs = make_variables(:tv, 1)
    fv = heat(uvs, pvs, tvs)
    js = FastDifferentiation.jacobian(fv, uvs) # The Jacobian js is of shape Ns x Ns

    jvps, vvs = FastDifferentiation.jacobian_transpose_v(fv, uvs) # The variable jvps is of shape Ns and rvs is of shape Ns+2

    input_vec = [uvs; pvs; tvs; vvs]
    jtv = make_function(FastDifferentiation.Node.(js'vvs), input_vec)
    jtv2 = make_function(jvps, input_vec)

    for i in 1:100
        rng = Random.Xoshiro(123)
        tvec = rand(rng, length(input_vec))
        @test isapprox(jtv(tvec), jtv2(tvec))
    end
end

@testitem "hessian" begin
    import FastDifferentiation as FD
    FD.@variables x y z

    h = hessian(x^2 * y^2 * z^2, [x, y, z])
    h_exe = FD.make_function(h, [x, y, z])

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
    import FastDifferentiation as FD

    FD.@variables x y z

    h = sparse_hessian(x^2 * y^2 * z^2, [x, y, z])
    h_exe = FD.make_function(h, [x, y, z])
    inp = sparse(ones(Float64, 3, 3))

    @test isapprox(
        h_exe([1, 2, 3]),
        [
            72.0 72.0 48.0
            72.0 18.0 24.0
            48.0 24.0 8.0]
    )

    h2_exe = FD.make_function(h, [x, y, z], in_place=true)
    h2_exe(inp, [1, 2, 3])
    @test isapprox(inp,
        [
            72.0 72.0 48.0
            72.0 18.0 24.0
            48.0 24.0 8.0])

    #test to make sure sparse_hessian/evaluate_path bug is not reintroduced (ref commit 4b4aeeb1990a15443ca87c15638dcaf7bd9d34d1)
    a = hessian(x * y, [x, y])
    b = sparse_hessian(x * y, [x, y])
    @test all(a .== b)
end

@testitem "hessian_times_v" begin
    using StaticArrays
    import FastDifferentiation as FD

    FD.@variables x y
    v_vec = make_variables(:v, 2)
    f = x^2 * y^2
    hv_slow = convert.(FD.Node, hessian(f, [x, y]) * v_vec)

    hv_slow_exe = FD.make_function(hv_slow, [[x, y]; v_vec])

    hv_fast, v_vec2 = hessian_times_v(f, [x, y])

    hv_fast_exe = FD.make_function(hv_fast, [[x, y]; v_vec2])

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
    import FastDifferentiation as FD

    FD.@variables a11 a12 a13 a21 a22 a23 a31 a32 a33

    vars = vec([a11 a12 a13 a21 a22 a23 a31 a32 a33])
    spmat = sparse([a11 a12 a13; a21 a22 a23; a31 a32 a33])
    f1 = FD.make_function(spmat, vars)
    inputs = [1 2 3 4 5 6 7 8 9]
    correct = [1 2 3; 4 5 6; 7 8 9]
    inp = similar(sprand(3, 3, 1.0))
    f2 = FD.make_function(spmat, vars, in_place=true)

    @test f1(inputs) == correct
    f2(inp, inputs)
    @test inp == correct
end

@testitem "SArray return" begin
    include("ShareTestCode.jl")
    using StaticArrays
    import FastDifferentiation as FD

    FD.@variables x y
    j = jacobian([x^2 * y^2, cos(x + y), log(x / y)], [x, y])
    j_exe = FD.make_function(j, [x, y])
    @test typeof(j_exe([1.0, 2.0])) <: Array
    j_exe2 = FD.make_function(SArray{Tuple{3,2}}(j), [x, y])
    @test typeof(j_exe2([1.0, 2.0])) <: StaticArray

    test_vec = [1.1, 2.3, 3.1]
    sph_func = FDTests.spherical_harmonics(4)
    sph_jac = jacobian(FD.roots(sph_func), FD.variables(sph_func))
    mn_func1 = FD.make_function(sph_jac, FD.variables(sph_func)) #return type of executable should be Array
    m, n = size(sph_jac)
    mn_func2 = FD.make_function(SMatrix{m,n}(sph_jac), FD.variables(sph_func)) #return type of executable should be StaticArray
    @test typeof(mn_func1(test_vec)) <: Array
    @test typeof(mn_func2(test_vec)) <: StaticArray

    @test isapprox(mn_func1(test_vec), mn_func2(test_vec))
    @test isapprox(mn_func1(SVector{3}(test_vec)), mn_func2(test_vec))
    @test isapprox(mn_func1(SVector{3}(test_vec)), mn_func2(SVector{3}(test_vec)))
    @test isapprox(mn_func1(test_vec), mn_func2(SVector{3}(test_vec)))
end

@testitem "dot bug and others" begin
    import FastDifferentiation as FD

    x = make_variables(:x, 2)
    mu = make_variables(:mu, 2)

    function no_exceptions()
        x'mu
        ex = sum(abs2, x .* mu)

        h = sparse_hessian(ex, x)
        hs = sparse_hessian(ex, x)


        fun = FD.make_function(h, [x; mu])
        fun = FD.make_function(hs, [x; mu])
        return true
    end

    @test no_exceptions()

    xprime = FD.Node.(x) #this will change the type of the vector
    fn = FD.make_function([xprime'mu], xprime, mu)

    @test isapprox(fn([1, 2, 3, 4])[1], 11)
end

@testitem "reverse_AD" begin
    include("ShareTestCode.jl")
    import FastDifferentiation as FD

    sph_func = FDTests.spherical_harmonics(7)
    sph_jac = jacobian(FD.roots(sph_func), FD.variables(sph_func))
    mn_func1 = FD.make_function(sph_jac, FD.variables(sph_func))

    rev_jac = similar(sph_jac)
    for (i, root) in pairs(FD.roots(sph_func))
        rev_jac[i, :] .= FD.reverse_AD(root, FD.variables(sph_func))
    end

    mn_func2 = FD.make_function(rev_jac, FD.variables(sph_func))

    test_vector = rand(3)
    @test isapprox(mn_func1(test_vector), mn_func2(test_vector), atol=1e-11)

end

@testitem "in place init_with_zeros" begin
    import FastDifferentiation as FD

    FD.@variables x y
    A = [x 0; 0 y]
    mat = [10 10; 10 10]

    fn = FD.make_function(A, [x, y], in_place=true, init_with_zeros=true)
    fn(mat, [1, 1])
    @test isapprox(mat, [1 0; 0 1])
    fn2 = FD.make_function(A, [x, y], in_place=true, init_with_zeros=false)
    mat = [10 10; 10 10]
    println(mat)
    fn2(mat, [1, 1])
    @test isapprox(mat, [1 10; 10 1])

    #NOT a test because of difficulty and fragility of parsing generated code. You have to verify these by looking at the output.
    p = make_variables(:p, 21)

    println("NO array zero statement")
    show(make_Expr(p, p, true, true))
    show(make_Expr(p, p, true, false))
    show(make_Expr(p, p, false, true))
    show(make_Expr(p, p, false, false))

    p[21] = 0

    println("shouldn't have an array zero statement but it should have a p[21]= 0 statement")
    show(make_Expr(p, p, true, true))
    println("this should not have an array zero statement nor should have a p[21] = 0 statement")
    show(make_Expr(p, p, true, false))
    println("should not have an array zero statement but should have a p[21] = 0 statement")
    show(make_Expr(p, p, false, true))
    show(make_Expr(p, p, false, false))

    p[20] = 0
    println("this should have an array zero statement should not have p[20]=0 or p[21]=0 statementt")
    show(make_Expr(p, p, true, true))
    println("this should not have an array zero statement should not have p[20]=0 or p[21]=0 statement")
    show(make_Expr(p, p, true, false))
    println("these should both have an array zero creation but should not have p[20]=0 or p[21]=0 statement")
    show(make_Expr(p, p, false, true))
    show(make_Expr(p, p, false, false))
end


@testitem "test c1*a ± c2*a => (c1 ± c1)*a" begin
    import FastDifferentiation as FD
    FD.@variables x

    f = 2x + 3x
    @test f === 5 * x
    f2 = -x + x
    @test FD.value(f2) == 0
    f3 = x + -x
    @test FD.value(f3) == 0
    f4 = 2x - 3x
    @test f4 === -x
    f5 = 2x - -3x
    @test f5 === 5x
end

@testitem "init with constants" begin
    using Random
    using SparseArrays
    using StaticArrays
    import FastDifferentiation as FD

    #difficult to test that inline array return code is being generated without unparsing generated function. But can at least verify that all different constant types of function arrays, all zeros, all non-zero constants, mix of the two, are filled properly.
    rng = Xoshiro(9398)
    size = 5
    for p in 0.0:0.1:1.0
        mat = convert(Matrix, sprand(rng, size, size, p))
        nmat2 = SMatrix{size,size,FD.Node}(FD.Node.(mat))
        nmat = FD.Node.(mat)

        println(nmat2)
        func1 = FD.make_function(nmat, FD.Node[])
        func2 = FD.make_function(nmat2, FD.Node[])
        println(func2)

        @assert isapprox(func1(Float64[]), mat)
        @assert isapprox(func2(Float64[]), mat)
    end
end

@testitem "more init with constant" begin
    import FastDifferentiation as FD

    function test_code_generation(f, input)
        correct_result = f(input)
        x_node = make_variables(:x, length(input))
        f_node = FD.Node.(f(x_node))

        # out of place
        f_callable = FD.make_function(f_node, x_node)
        @show f_callable
        result = f_callable(input)
        @test isapprox(result, correct_result)
        @test typeof(result) <: typeof(correct_result)

        # in_place, init_with_zeros
        f_callable_init! = FD.make_function(f_node, x_node; in_place=true, init_with_zeros=true)
        @show f_callable_init!
        result = similar(correct_result)
        f_callable_init!(result, input) == correct_result
        @test isapprox(result, correct_result)

        # in_place, !init_with_zeros
        f_callable_no_init! = FD.make_function(f_node, x_node; in_place=true, init_with_zeros=false)
        @show f_callable_no_init!
        result = similar(correct_result)
        result_copy = copy(result)
        f_callable_no_init!(result, input)
        for i in eachindex(result)
            if !iszero(correct_result[i])
                @test isapprox(result[i], correct_result[i])
            end
        end
    end

    # systematically enumerate the four different branches of `make_Expr`:
    # all constant, mostly zeros
    @info "all constant, mostly zeros"
    test_code_generation([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]) do x
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    end

    # all constant, some zeros
    @info "all constant, some zeros"
    test_code_generation([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]) do x
        [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0]
    end

    # mostly constants
    @info "mostly constants"
    test_code_generation(3.0) do x
        [2.1 * x[1], 1, 2]
    end

    # non-constant at non-first position
    @info "non-constant at non-first position"
    test_code_generation(3.0) do x
        [1, 2.1 * x[1], 2]
    end

    # mostly zeros
    @info "mostly zeros"
    test_code_generation(3.0) do x
        [x[1]^2, 0, 0]
    end

    # all non-constant
    @info "all non-constant"
    test_code_generation(3.0) do x
        [2.1 * x[1], x[1]^2, sqrt(x[1])]
    end

    @info "evaluate with exotic eltype"
    test_code_generation(Complex(1.0)) do x
        [1, 2.1 * x[1], 2]
    end
end

@testitem "no variables in function" begin
    @variables x y

    f = x + y
    h = hessian(f, [x, y])
    hexe = make_function(h, [x, y])

    @test isapprox(hexe([1, 1]), [0 0; 0 0])
end



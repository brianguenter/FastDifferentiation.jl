module FDTests
using StaticArrays
using Memoize
using DataStructures

using FastDifferentiation



const FD = FastDifferentiation
export FD

include("TestPrograms/TestCode.jl")

"""If `compute_dominators` is `true` then computes `idoms` tables for graph, otherwise computes `pidoms` table`"""
function compute_dominance_tables(graph::FD.DerivativeGraph{T}, compute_dominators::Bool) where {T<:Integer}
    if compute_dominators
        start_vertices = FD.root_index_to_postorder_number(graph)
    else
        start_vertices = FD.variable_index_to_postorder_number(graph)
    end

    doms = Dict{T,T}[]   #create one idom table for each root

    for (start_index, node_postorder_number) in pairs(start_vertices)
        push!(doms, FastDifferentiation.compute_dom_table(graph, compute_dominators, start_index, node_postorder_number))
    end
    return doms
end
export compute_dominance_tables


function simple_dag(cache::Union{IdDict,Nothing}=IdDict())
    FD.@variables zz y
    return expr_to_dag(zz^2 + y * (zz^2), cache), zz, y
end
export simple_dag

function simple_numbered_dag()
    FD.@variables zz
    return expr_to_dag(zz * cos(zz^2))
end
export simple_numbered_dag

function dominators_dag()
    FD.@variables zz
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
    FD.@variables x

    sinx = FD.Node(sin, MVector(x))
    cosx = FD.Node(cos, MVector(x))
    A = FD.Node(*, MVector(cosx, sinx))
    sinA = FD.Node(sin, MVector(A))
    expsin = FD.Node(*, MVector(A, FD.Node(exp, MVector(sinx))))
    plus = FD.Node(+, MVector(sinA, expsin))
    return plus
end
export complex_dominator_dag

complex_dominator_graph() = FD.DerivativeGraph(complex_dominator_dag())
export complex_dominator_graph

function R2_R2_function()
    FD.@variables x y

    n2 = x * y
    n4 = n2 * y
    n5 = n2 * n4

    return FD.DerivativeGraph([n5, n4])
end
export R2_R2_function

function simple_dominator_graph()
    FD.@variables x

    ncos = FD.Node(cos, x)
    nplus = FD.Node(+, ncos, x)
    ntimes = FD.Node(*, ncos, nplus)
    four_2_subgraph = FD.Node(+, nplus, ncos)
    one_3_subgraph = FD.Node(+, FD.Node(*, FD.Node(-1), FD.Node(sin, x)), FD.Node(1))
    return x, FD.DerivativeGraph(ntimes), four_2_subgraph, one_3_subgraph
end
export simple_dominator_graph

"""returns 4 factorable subgraphs in this order: (4,2),(1,3),(1,4),(4,1)"""
function simple_factorable_subgraphs()
    _, graph, _, _ = simple_dominator_graph()
    temp = extract_all!(FD.compute_factorable_subgraphs(graph))
    return graph, [
        temp[findfirst(x -> FastDifferentiation.vertices(x) == (4, 2), temp)],
        temp[findfirst(x -> FastDifferentiation.vertices(x) == (1, 3), temp)],
        temp[findfirst(x -> FastDifferentiation.vertices(x) == (1, 4), temp)],
        temp[findfirst(x -> FastDifferentiation.vertices(x) == (4, 1), temp)]
    ]
end
export simple_factorable_subgraphs
end #module

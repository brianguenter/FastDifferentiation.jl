module Vis

using Colors
using ColorTypes
using FixedPointNumbers
using ..FastSymbolicDifferentiation

using GraphRecipes
using Plots
using Graphs

# digraph graphname {
#     // This attribute applies to the graph itself
#     size="1,1";
#     // The label attribute can be used to change the label of a node
#     a [label="Foo"];
#     // Here, the node shape is changed.
#     b [shape=box];
#     // These edges both have different line properties
#     a -> b -> c [color=blue];
#     b -> d [style=dotted];
#     // [style=invis] hides a node.
# }
function to_dot(graph::RnToRmGraph, filename)
    data = """
    graph graphname {
        // This attribute applies to the graph itself
         size="1,1";
         // The label attribute can be used to change the label of a node
         a [label="Foo"];
         // Here, the node shape is changed.
         b [shape=box];
         // These edges both have different line properties
         a -- b -- c [color=blue];
         b -- d [style=dotted];
         // [style=invis] hides a node.
    }"""
    f = open(filename, "w")
    write(f, data)
    close(f)
end
export to_dot

function explot()
    n = 8
    g = wheel_digraph(n)
    edgelabel_dict = Dict()
    edgelabel_mat = Array{String}(undef, n, n)
    for j in 1:n
        for i in 1:n
            # edgelabel_mat[i, j] = 
            edgelabel_dict[(i, j)] = string("edge ", i, " to ", j)
        end
    end
    # edgelabel_vec = edgelabel_mat[:]

    graphplot(g, names=1:n, edgelabel=edgelabel_dict, curves=false, nodeshape=:rect)  # Or edgelabel=edgelabel_mat, or edgelabel=edgelabel_vec.
end
export explot

function graph_recipes_example()
    # default(size=(2000, 2000))

    #         out += 1
    #     end
    #     out
    # end
    # )

    code = dag_to_Expr(complex_dominator_dag(), Dict{Node,Union{Symbol,Real}}())

    plot(code, fontsize=7, shorten=0.01, axis_buffer=0.15, nodeshape=:rect, nodecolor="pink")
end
export graph_recipes_example

function label_func(mask::BitVector, label_string::String)
    root_string = ""
    for (i, bvar) in pairs(mask)
        if bvar == 1
            root_string *= "$label_string$i,"
        end
    end
    return root_string
end

function edge_label(a::PathEdge)
    root_string = label_func(reachable_roots(a), "r")
    var_string = label_func(reachable_variables(a), "v")
    return ("$root_string : $var_string")
end
export edge_label

function draw(unique_edges::AbstractVector{PathEdge{T}}, unique_nodes::AbstractVector{Node{T,N}}, post_order_numbers::AbstractVector{T}, root_test::Function) where {T,N}
end

"""draws nodes and labled edges of a DerivativeGraph"""
function draw(graph, value_labels=true; draw_edge_labels=true, draw_node_labels=true)
    default(size=(1000, 1200))
    unique_nodes = nodes(graph)
    num_nodes = length(unique_nodes)
    g = SimpleGraph(num_nodes)
    node_labels = fill("", num_nodes) #Vector{String}(undef, num_nodes)
    edge_labels = Dict()
    node_fill = RGBA{Float64}[]
    green = RGBA(0.0, 1.0, 0.0)
    purple = RGBA(1.0, 0.0, 0.0)
    lightblue = RGBA(0.4, 0.4, 0.0)
    yellow = RGBA(1.0, 1.0, 0.0, 0.0)


    for gnode in unique_nodes
        node_index = postorder_number(graph, gnode)
        if is_root(graph, node_index)
            push!(node_fill, green)
        elseif is_variable(graph, node_index)
            push!(node_fill, purple)
        elseif is_constant(graph, node_index)
            push!(node_fill, yellow)
        else
            push!(node_fill, lightblue)
        end

        if draw_node_labels
            node_labels[node_index] = "$(postorder_number(graph,gnode)): $(node_value(gnode))"
        end
    end

    un = unique_edges(graph)
    tmp_edges = collect(un)
    multi_edges = Vector{Any}[]

    #find multiple edges between nodes
    for edge in tmp_edges
        sub = (top_vertex(edge), bott_vertex(edge))
        medges = filter(x -> (top_vertex(x), bott_vertex(x)) == sub, un)

        if length(medges) != 0
            delete!.(Ref(un), medges)
            push!(multi_edges, collect(medges))
        end
    end

    for edge_group in multi_edges
        first_edge = edge_group[1]
        add_edge!(g, bott_vertex(first_edge), top_vertex(first_edge))

        if draw_edge_labels
            label = ""

            for edge in edge_group
                if !value_labels
                    label *= "(" * edge_label(edge) * ") -- "
                else
                    label *= "(" * (node_value(edge_value(edge)) isa AutomaticDifferentiation.NoDeriv ? "NoDeriv" : string(edge_value(edge))) * ") -- "
                end
            end
            edge_labels[(bott_vertex(first_edge), top_vertex(first_edge))] = string(label)
        else
            edge_labels = String[]
        end
    end

    edge_colors = [RGBA(0.0, 0.0, 0.0, 0.0) for i in 1:10]
    # edge_colors = [RGBA(1.0, 0.5, 0.0, 0.7) for i in 1:10]
    tmp = graphplot(g,
        curves=false,
        edgelabel=edge_labels,
        names=node_labels,
        nodelabelc="white",
        edgelabelc="green",
        # markercolor="lightblue",
        nodecolor=node_fill,
        nodeshape=:rect,
        edgestrokec=edge_colors,
        linecolor=:red,
        method=:tree,
        fontsize=10,
        nodesize=0.05,
        markerstrokewidth=0,
        edgelabelbox=true,
    )
    display(tmp)
    # readline()
end
export draw
end #module
export Vis

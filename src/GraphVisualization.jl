module Vis

using Colors
using ColorTypes
using ..FastSymbolicDifferentiation

using GraphRecipes
using Plots
using Graphs

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

# function draw(unique_edges::AbstractVector{PathEdge{T}}, unique_nodes::AbstractVector{Node{T,N}}, post_order_numbers::AbstractVector{T}, root_test::Function) where {T,N}
# end

function draw_dot(graph, filename)
    gr = "strict digraph{\nnode [style = filled]\n"
    for e in FastSymbolicDifferentiation.unique_edges(graph)
        roots = join(findall(x -> x == 1, reachable_roots(e)), ",")
        variables = join(findall(x -> x == 1, reachable_variables(e)), ",")
        gr *= "$(top_vertex(e)) -> $(bott_vertex(e)) [label = \"r:[$roots]  v:[$variables]\"] [color = purple]\n"
    end

    for node in nodes(graph)
        num = postorder_number(graph, node)
        if is_variable(graph, num)
            gr *= "$num [color = green] [label = \"v$(variable_postorder_to_index(graph,num)) $(value(node))\"] [fillcolor = \"#96ff96\"]\n"
        elseif is_root(graph, num)
            gr *= "$num [color = red] [label = \"r$(root_postorder_to_index(graph,num)) $num $(value(node))\"] [fillcolor = \"#ff9696\"]\n"
        elseif is_constant(graph, num)
            gr *= "$num [color = \"#969600\"] [label = \"$(value(node))\"] [fillcolor = \"#ffff00\"]\n"
        else
            gr *= "$num [label = \"$num $(value(node))\"]\n"
        end
    end

    gr *= "\n}"
    name, ext = splitext(filename)
    io = open(name * ".dot", "w")
    write(io, gr)
    close(io)

    Base.run(`dot -Tsvg $name.dot -o $filename`)
end
export draw_dot


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
            node_labels[node_index] = "$(postorder_number(graph,gnode)): $(value(gnode))"
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
                    label *= "(" * (value(value(edge)) isa AutomaticDifferentiation.NoDeriv ? "NoDeriv" : string(value(edge))) * ") -- "
                end
            end
            edge_labels[(bott_vertex(first_edge), top_vertex(first_edge))] = string(label)
        else
            edge_labels = String[]
        end
    end

    println(sum(length.(g.fadjlist)))

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
        fontsize=4,
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

module Vis

using Colors
using ColorTypes
using ..FastSymbolicDifferentiation
using ..FastSymbolicDifferentiation: PathEdge, nodes, postorder_number, is_root, is_variable, is_constant, value, unique_edges, top_vertex, bott_vertex, AutomaticDifferentiation, reachable_roots, reachable_variables, node, parent_edges, variable_postorder_to_index, root_postorder_to_index, DerivativeGraph
using GraphRecipes
using Plots
using Graphs
using ElectronDisplay

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


# function draw(unique_edges::AbstractVector{PathEdge{T}}, unique_nodes::AbstractVector{Node{T,N}}, post_order_numbers::AbstractVector{T}, root_test::Function) where {T,N}
# end

function edges_from_node(graph, start_node_number::AbstractVector{Int})
    result = Set{PathEdge}()
    work_queue = Int[]
    append!(work_queue, start_node_number)
    while length(work_queue) != 0
        curr_node = pop!(work_queue)
        new_edges = child_edges(graph, curr_node)
        if length(new_edges) != 0
            push!(result, new_edges...)
        end
        for edge in new_edges
            if !in(bott_vertex(edge), work_queue)
                push!(work_queue, bott_vertex(edge))
            end
        end
    end
    return collect(result)
end

function make_dot_file(graph, start_nodes::Union{Nothing,AbstractVector{Int}}, label::String, reachability_labels=true, value_labels=true, no_path_edges=false)
    if start_nodes !== nothing
        edges_to_draw = edges_from_node(graph, start_nodes)
    else
        edges_to_draw = collect(FastSymbolicDifferentiation.unique_edges(graph))
    end

    return make_dot_file(graph, edges_to_draw, label, reachability_labels, value_labels, no_path_edges)

end
export make_dot_file

function make_dot_file(graph, edges_to_draw::AbstractVector{P}, label::String, reachability_labels=true, value_labels=true, no_path_edges=false) where {P<:PathEdge}
    if !no_path_edges #only draw edges on path from root to variable
        edges_to_draw = collect(filter(x -> any(reachable_variables(x)) && any(reachable_roots(x)), edges_to_draw))
    end

    gr = "digraph{\nnode [style = filled]\n"
    if label != ""
        gr *= "label = \"$label\"\n"
    end
    gr *= "ratio=\"fill\"\n"
    # gr *= "size = 12 12\n"
    gr_copy = deepcopy(graph)

    nodes_to_draw = Set{Node}()
    for e in edges_to_draw
        roots = join(findall(x -> x == 1, reachable_roots(e)), ",")
        variables = join(findall(x -> x == 1, reachable_variables(e)), ",")
        edge_label = ""

        if value_labels
            edge_label *= "$(value(e)) "
        end
        if reachability_labels
            edge_label *= "  r:[$roots]  v:[$variables]"
        end

        if edge_label != ""
            edge_label = "[label = " * "\"" * edge_label * "\"] [color = purple]"
        end


        gr *= "$(top_vertex(e)) -> $(bott_vertex(e)) $edge_label\n"

        push!(nodes_to_draw, node(gr_copy, top_vertex(e)))
        push!(nodes_to_draw, node(gr_copy, bott_vertex(e)))
    end

    for node in nodes_to_draw
        if !(!is_root(gr_copy, node) && length(parent_edges(gr_copy, node)) == 0 && length(child_edges(gr_copy, node)) == 0)
            num = postorder_number(gr_copy, node)
            if is_variable(gr_copy, num)
                gr *= "$num [color = green] [label = \"v$(variable_postorder_to_index(gr_copy,num)) $num $(value(node))\"] [fillcolor = \"#96ff96\"]\n"
            elseif is_root(gr_copy, num)
                gr *= "$num [color = red] [label = \"r$(root_postorder_to_index(gr_copy,num)) $num $(value(node))\"] [fillcolor = \"#ff9696\"]\n"
            elseif is_constant(gr_copy, num)
                gr *= "$num [color = \"#969600\"] [label = \"$num  $(value(node))\"] [fillcolor = \"#ffff00\"]\n"
            else
                gr *= "$num [label = \"$num $(value(node))\"]\n"
            end
        end
    end

    gr *= "\n}"
    return gr
end

function draw_dot(graph, edges_to_draw::AbstractVector{P}; graph_label::String="", reachability_labels=true, value_labels=true) where {P}
    gr = make_dot_file(graph, edges_to_draw, graph_label, reachability_labels, value_labels)
    path, io = mktemp(cleanup=true)
    name, ext = splitext(path)
    write_dot(name * ".svg", gr)
    svg = read(name * ".svg", String)
    display("image/svg+xml", svg)
end
export draw_dot

function draw_dot(graph; start_nodes::Union{Nothing,AbstractVector{Int}}=nothing, graph_label::String="", reachability_labels=true, value_labels=true)
    path, io = mktemp(cleanup=true)
    name, ext = splitext(path)
    write_dot(name * ".svg", graph;
        start_nodes=start_nodes,
        graph_label=graph_label,
        reachability_labels=reachability_labels,
        value_labels=value_labels,
        no_path_edges=true)
    svg = read(name * ".svg", String)
    display("image/svg+xml", svg)
end
export draw_dot

draw_dot(subgraph::FastSymbolicDifferentiation.FactorableSubgraph, graph_label::String="", reachability_labels=true, value_labels=false) = draw_dot(graph(subgraph), collect(subgraph_edges(subgraph)), graph_label=graph_label, reachability_labels=reachability_labels, value_labels=value_labels)

function write_dot(filename, graph::DerivativeGraph; start_nodes::Union{Nothing,AbstractVector{Int}}=nothing, graph_label::String="", reachability_labels=true, value_labels=true, no_path_edges=false)
    gr = make_dot_file(graph, start_nodes, graph_label, reachability_labels, value_labels, no_path_edges)
    write_dot(filename, gr)
end
export write_dot

function write_dot(filename, graph::DerivativeGraph, edges_to_draw::AbstractVector{PathEdge}; start_nodes::Union{Nothing,AbstractVector{Int}}=nothing, graph_label::String="", reachability_labels=true, value_labels=true, no_path_edges=false)
    gr = make_dot_file(graph, edges_to_draw, graph_label, reachability_labels, value_labels, no_path_edges)
    write_dot(filename, gr)
end
export write_dot
function write_dot(filename, dot_string::String)
    name, ext = splitext(filename)
    name = name * ".dot"
    ext = ext[2:end] #get rid of leading dot
    io = open(name, "w")
    write(io, dot_string)
    close(io)

    Base.run(`dot -T$ext $name -o $filename`)
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

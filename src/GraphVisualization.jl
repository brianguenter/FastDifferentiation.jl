module Vis

using Colors
using ColorTypes
using ..FastSymbolicDifferentiation

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
export edge_label

# function draw(unique_edges::AbstractVector{PathEdge{T}}, unique_nodes::AbstractVector{Node{T,N}}, post_order_numbers::AbstractVector{T}, root_test::Function) where {T,N}
# end

function edges_from_node(graph, start_node_number::AbstractVector{Int})
    result = PathEdge[]
    work_queue = Int[]
    append!(work_queue, start_node_number)
    while length(work_queue) != 0
        curr_node = pop!(work_queue)
        new_edges = child_edges(graph, curr_node)
        append!(result, new_edges)
        for edge in new_edges
            if !in(bott_vertex(edge), work_queue)
                push!(work_queue, bott_vertex(edge))
            end
        end
    end
    return result
end

function make_dot_file(graph, start_nodes::Union{Nothing,AbstractVector{Int}}, label::String, reachability_labels=true, value_labels=true)
    gr = "strict digraph{\nnode [style = filled]\n"
    if label != ""
        gr *= "label = \"$label\"\n"
    end
    gr *= "ratio=\"fill\"\n"
    gr *= "size = 12 12\n"
    gr_copy = deepcopy(graph)
    FastSymbolicDifferentiation.remove_dangling_edges!(gr_copy)

    if start_nodes !== nothing
        edges_to_draw = edges_from_node(gr_copy, start_nodes)
    else
        edges_to_draw = FastSymbolicDifferentiation.unique_edges(gr_copy)
    end

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
export make_dot_file

function draw_dot(graph; start_nodes::Union{Nothing,AbstractVector{Int}}=nothing, graph_label::String="", reachability_labels=true)
    gr = make_dot_file(graph, start_nodes, graph_label, reachability_labels)
    path, io = mktemp(cleanup=true)
    name, ext = splitext(path)
    write(io, gr)
    close(io)
    println(name)
    Base.run(`dot -Tsvg $path -o $name.svg`)
    svg = read(name * ".svg", String)
    display("image/svg+xml", svg)
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

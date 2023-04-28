# module GraphProcessing
# using FastSymbolicDifferentiation: RnToRmGraph, roots, variables, nodes, _node_edges, children, parents, root_index_to_postorder_number, variable_index_to_postorder_number, each_vertex, top_vertex, bott_vertex, PathEdge, edges, reachable_roots, reachable_variables, codomain_dimension, domain_dimension, subset, is_zero, to_string, EdgeRelations
using StaticArrays
using DataStructures


"""This and PathConstraint are highly redundant. There is one function, compute_dominance_tables, which used PathConstraint in a slightly different way than all the other code. This caused endless bugs. For now make a new constraint struct just for compute_dominance_tables, since that code is complicated and don't want to rewrite it now. Later get rid of DomPathConstraint struct, since it is only used by one function."""
struct DomPathConstraint{T<:Integer}
    graph::DerivativeGraph{T}
    iterate_parents::Bool
    roots_mask::BitVector
    variables_mask::BitVector
    relations::Vector{T} #scratch vector that is used to hold temporary values for iteration


    function DomPathConstraint(graph::DerivativeGraph{T}, iterate_parents::Bool, roots_mask::BitVector, variables_mask::BitVector) where {T}
        relations = Vector{T}(undef, 5) #initialize to small size since num parents or children is likely to be small for any given node. 
        # @assert !is_zero(roots_mask) #can't iterate with all roots constrained out
        # @assert !is_zero(variables_mask) #can't iterate with all variables constrained out
        return new{T}(graph, iterate_parents, roots_mask, variables_mask, relations)
    end
end
export DomPathConstraint

function DomPathConstraint(graph::DerivativeGraph, iterate_parents::Bool, root_or_leaf_index::Integer)
    if iterate_parents
        roots_mask = falses(codomain_dimension(graph))
        roots_mask[root_or_leaf_index] = 1
        variables_mask = falses(domain_dimension(graph))
    else
        roots_mask = falses(codomain_dimension(graph))
        variables_mask = falses(domain_dimension(graph))
        variables_mask[root_or_leaf_index] = 1
    end

    return DomPathConstraint(graph, iterate_parents, roots_mask, variables_mask)
end

graph(a::DomPathConstraint) = a.graph
roots_mask(a::DomPathConstraint) = a.roots_mask
variables_mask(a::DomPathConstraint) = a.variables_mask
relations(a::DomPathConstraint) = a.relations

"""Contains information used to constrain graph traversal to only edges/vertices that are on the path to a root or variable vertex"""
struct PathConstraint{T<:Integer}
    dominating_node::T
    graph::DerivativeGraph{T}
    iterate_parents::Bool
    roots_mask::BitVector
    variables_mask::BitVector
    relations::Vector{T} #scratch vector that is used to hold temporary values for iteration


    function PathConstraint(dominating_node::T, graph::DerivativeGraph{T}, iterate_parents::Bool, roots_mask::BitVector, variables_mask::BitVector) where {T}
        relations = Vector{T}(undef, 5) #initialize to small size since num parents or children is likely to be small for any given node. 
        # @assert !is_zero(roots_mask) #can't iterate with all roots constrained out
        # @assert !is_zero(variables_mask) #can't iterate with all variables constrained out
        return new{T}(dominating_node::T, graph, iterate_parents, roots_mask, variables_mask, relations)
    end
end
export PathConstraint

function PathConstraint(graph::DerivativeGraph, iterate_parents::Bool, root_or_leaf_index::Integer)
    if iterate_parents
        roots_mask = falses(codomain_dimension(graph))
        roots_mask[root_or_leaf_index] = 1
        variables_mask = falses(domain_dimension(graph))
    else
        roots_mask = falses(codomain_dimension(graph))
        variables_mask = falses(domain_dimension(graph))
        variables_mask[root_or_leaf_index] = 1
    end

    return PathConstraint(graph, iterate_parents, roots_mask, variables_mask)
end

dominating_node(a::PathConstraint) = a.dominating_node
roots_mask(a::PathConstraint) = a.roots_mask
variables_mask(a::PathConstraint) = a.variables_mask

# path_mask(a::PathConstraint) = a.roots_mask
# export path_mask
graph(a::PathConstraint) = a.graph
export graph


is_dominator_constraint(a::PathConstraint) = a.iterate_parents

function Base.show(io::IO, a::PathConstraint)
    print(io, "iterate_parents $(a.iterate_parents) [$(to_string(roots_mask(a),"r")) ↔ $(to_string(variables_mask(a),"v"))]")
end

"""Returns node indices connected to `node_index` which satisfy the path constraint. This is different from edge_relations functions below which return edges, not node indices."""
function relation_node_indices(a::DomPathConstraint{T}, node_index::T) where {T<:Integer}
    # node_relations = a.relations #use scratch array in path constraint. This should only be used by a single thread so this should be multi-thread safe and avoids allocating many small arrays.
    # #TODO might want to change this. This seems intrinsically dangerous since node_relations is returned but is also a component of a PathConstraint. Could lead to very subtle bugs. The only caller of this function is compute_dominance_tables. The path constraint is created in a loop and used and consumed in that loop so it should never be possible to another piece of code to mess with this variable. Still seems too tricky.
    # empty!(node_relations) #reset array to zero
    node_relations = T[] #for now do not reuse a.relations since this could cause subtle bugs. Optimize later

    curr_edges = node_edges(graph(a), node_index)
    if curr_edges === nothing
        return nothing
    end

    #these tests allow for either empty variables_mask or roots_mask. The functions which call this depend on being able to do this.
    if a.iterate_parents
        for edge in parents(curr_edges)
            if (bott_vertex(edge) == node_index) && subset(roots_mask(a), reachable_roots(edge)) && subset(variables_mask(a), reachable_variables(edge)) #top_vertex(edge) is a parent of node_index and top_vertex(edge) is on the path to the correct root or leaf
                push!(node_relations, top_vertex(edge))
            end
        end
    else
        for edge in children(curr_edges)
            if (top_vertex(edge) == node_index) && subset(variables_mask(a), reachable_variables(edge)) && subset(roots_mask(a), reachable_roots(edge)) #bott_vertex(edge) is a child of node_index and bott_vertex(edge) is on the path to the correct root or leaf
                push!(node_relations, bott_vertex(edge))
            end
        end
    end

    if length(node_relations) == 0
        return nothing
    else
        return node_relations #this is kind of bad because functions can mess with relations instance variable. But ConstrainedPathIterator doesn't need it to remain unchanged between calls to relations. 
    end
end
export relation_node_indices


# struct RelationIterator{T<:Integer}
#     constraint::PathConstraint{T}
#     edges_of_node::EdgeRelations{T}
# end

# function iterate(a::RelationIterator) 
#     tmp_edges = _node_edges(graph_edges(a), node_index)
#     if a.iterate_parents
#         for edge in parents(tmp_edges)
#             if subset(roots_mask(a), reachable_roots(edge)) && any(variables_mask(a) .& reachable_variables(edge))
#                 push!(result, edge)
#             end
#         end
#     else
#         for edge in children(tmp_edges)
#             if subset(variables_mask(a), reachable_variables(edge)) && any(roots_mask(a) .& reachable_roots(edge))
#                 push!(result, edge)
#             end
#         end
#     end
# end

"""returns edges emanating from a vertex which satisfy the PathConstraint"""
function relation_edges!(a::PathConstraint{T}, node_index::Integer, result::Union{Nothing,Vector{PathEdge{Int64}}}=nothing) where {T<:Integer}
    if result === nothing
        result = PathEdge{T}[]
    else
        empty!(result)
    end

    tmp_edges = node_edges(graph(a), node_index)
    if tmp_edges === nothing
        return nothing
    end

    #these tests allow for either empty variables_mask or roots_mask. The functions which call this depend on being able to do this.
    if a.iterate_parents
        for edge in parents(tmp_edges)
            if subset(roots_mask(a), reachable_roots(edge)) && overlap(variables_mask(a), reachable_variables(edge))
                push!(result, edge)
            end
        end
    else
        for edge in children(tmp_edges)
            if subset(variables_mask(a), reachable_variables(edge)) && overlap(roots_mask(a), reachable_roots(edge))
                push!(result, edge)
            end
        end
    end

    return result
end

"""returns edges along a single path which satisfy the PathConstraint"""
function relation_edges!(a::PathConstraint{T}, edge::PathEdge, result::Union{Nothing,Vector{PathEdge{Int64}}}=nothing) where {T<:Integer}
    if result === nothing
        result = PathEdge{T}[]
    else
        empty!(result)
    end
    #these tests do not allow for empty variables_mask or roots_mask. The functions which call this variant of this function depend on this. Hacky, should be fixed. Later.
    if a.iterate_parents
        tmp = node_edges(graph(a), top_vertex(edge))

        if tmp === nothing
            return nothing
        end

        tmp_edges = parents(tmp)
        for one_edge in tmp_edges
            if top_vertex(one_edge) ≤ dominating_node(a) && subset(roots_mask(a), reachable_roots(one_edge)) && overlap(variables_mask(a), reachable_variables(one_edge))
                push!(result, one_edge)
            end
        end
    else
        tmp = node_edges(graph(a), bott_vertex(edge))

        if tmp === nothing
            return nothing
        end

        tmp_edges = children(tmp)

        for one_edge in tmp_edges
            if bott_vertex(one_edge) ≥ dominating_node(a) && subset(variables_mask(a), reachable_variables(one_edge)) && overlap(roots_mask(a), reachable_roots(one_edge))
                push!(result, one_edge)
            end
        end
    end
    return result
end
export relation_edges!

function _compute_paths!(path_masks, graph_edges::Dict{T,EdgeRelations{T}}, current_node_index, origin_index, relation_function) where {T}
    @assert current_node_index <= length(path_masks)
    if path_masks[current_node_index][origin_index] == 1
        return #already visited this node so don't recurse.
    else
        path_masks[current_node_index][origin_index] = 1
        nodes = relation_function(graph_edges, current_node_index)
        if nodes !== nothing #can have a constant or a variable as a root node so there won't be any children. Don't follow paths in this case.
            for next_node in nodes
                _compute_paths!(path_masks, graph_edges, next_node, origin_index, relation_function)
            end
        end
    end
end


"""Returns a vector of BitArray. The BitArray encodes the reachability of root (or leaf) nodes from a node `nᵢ`. Call the vector of BitArray `path_masks`: then `path_masks[i]` is the path mask for node `nᵢ`. If root(or leaf) `j` is reachable from node `nᵢ` then `path_masks[i][j] = 1` , 0 otherwise."""
function compute_paths(num_nodes::Integer, graph_edges::Dict{T,EdgeRelations{T}}, origin_nodes::Vector{T}, relation_function::Function) where {T<:Integer}
    path_masks = [falses(length(origin_nodes)) for _ in 1:num_nodes]

    for ((origin_index, postorder_num)) in pairs(origin_nodes)
        @assert postorder_num <= num_nodes
        _compute_paths!(path_masks, graph_edges, postorder_num, origin_index, relation_function)
    end

    return path_masks  #WARNING:TODO rewrite this code so don't create the array of path_masks. It is discarded after returning from the calling function compute_edge_paths! so this is wasteful and probably slow.
end

compute_paths_to_variables(num_nodes::Integer, graph_edges::Dict{T,EdgeRelations{T}}, var_indices) where {T} = compute_paths(num_nodes, graph_edges, var_indices, parents) #TODO: don't want to be converting to an array of integers everytime.
export compute_paths_to_variables
compute_paths_to_roots(num_nodes::Integer, graph_edges::Dict{T,EdgeRelations{T}}, indices_of_roots) where {T} = compute_paths(num_nodes, graph_edges, indices_of_roots, children)
export compute_paths_to_roots

"""This function computes reachable roots and reachable variable masks for each edge in the graph. It is called before the DerivativeGraph has been constructed. If you have a DerivativeGraph you can use compute_edge_paths!(a::DerivativeGraph) instead"""
function compute_edge_paths!(num_nodes::Integer, graph_edges::Dict{T,EdgeRelations{S}}, var_indices, indices_of_roots) where {T,S<:Integer}
    variable_paths = compute_paths_to_variables(num_nodes, graph_edges, var_indices)
    root_paths = compute_paths_to_roots(num_nodes, graph_edges, indices_of_roots)

    for node_index in keys(graph_edges)
        tmp_edges = _node_edges(graph_edges, node_index)
        for edge in children(tmp_edges) #this does twice as much work and allocates twice as much memory as it needs to because each edge appears twice in the edges of graph. Optimize later if necessary.
            copy!(edge.reachable_roots, root_paths[node_index])
        end
        for edge in parents(tmp_edges)
            copy!(edge.reachable_variables, variable_paths[node_index])
        end
    end
end
export compute_edge_paths!

"""convenience function to avoid calling messier low level function"""
compute_edge_paths!(graph::DerivativeGraph) = compute_edge_paths!(length(nodes(graph)), edges(graph), variable_index_to_postorder_number(graph), root_index_to_postorder_number(graph))


function intersection(order_test, node1::Integer, node2::Integer, idoms::Dict{T,T}) where {T<:Integer}
    count = 0
    max_count = length(idoms)
    while true
        #added this assertion because several simple errors in other code made this code loop forever. Hard to track the error down without a thrown exception.
        @assert count <= max_count "intersection has taken more steps than necessary. This should never happen."
        if node1 == node2
            return node1
        else
            if order_test(node1, node2)
                node1 = idoms[node1]
            else
                node2 = idoms[node2]
            end
        end
        count += 1
    end
end

function fill_idom_table!(next_vertices::Union{Nothing,AbstractVector{T}}, dom_table::Dict{T,T}, current_node::T, order_test::Function) where {T<:Integer}
    if next_vertices === nothing
        dom_table[current_node] = current_node
    elseif length(next_vertices) == 1
        dom_table[current_node] = next_vertices[1]
    else
        br1 = next_vertices[1]
        for relation_vertex in view(next_vertices, 2:length(next_vertices))
            br1 = intersection(order_test, br1, relation_vertex, dom_table)
        end
        dom_table[current_node] = br1
    end
end

"""compute the idom or pidom table starting from root or variable corresponding to start_index."""
function compute_dom_table(graph::DerivativeGraph{T}, compute_dominators::Bool, start_index::T, node_postorder_number::T, current_dom::Union{Nothing,Dict{T,T}}=nothing) where {T}
    if compute_dominators
        next_vertices_relation = (gr, curr_node::Integer) -> children(gr, curr_node)
        upward_path = true
        order_test = <
    else
        next_vertices_relation = (gr, curr_node::Integer) -> parents(gr, curr_node)
        upward_path = false
        order_test = >
    end

    path_constraint = DomPathConstraint(graph, compute_dominators, start_index)

    if current_dom === nothing
        current_dom = Dict{T,T}()
    else
        empty!(current_dom)
    end

    #this is only necessary if trying to multithread, otherwise can have single work_heap allocated outside for loop which is almost certainly more efficient for single threaded code.
    if compute_dominators
        work_heap = MutableBinaryMaxHeap{T}()
    else
        work_heap = MutableBinaryMinHeap{T}()
    end

    visited = Set{T}()

    push!(work_heap, node_postorder_number)

    #do BFS traversal of graph from largest postorder numbers downward. Don't think BFS is necessary and would probably be faster to use DFS without work_heap.


    while length(work_heap) != 0
        curr_level = length(work_heap)

        for _ in 1:curr_level

            curr_node = pop!(work_heap)
            @timeit TIMER "relation_node_indices" parent_vertices = relation_node_indices(path_constraint, curr_node) #for dominator this will return the parents of the current node, constrained to lie on the path to the start_vertex.

            @timeit TIMER "fill_idom_table" fill_idom_table!(parent_vertices, current_dom, curr_node, order_test)

            if next_vertices_relation(graph, curr_node) !== nothing
                #get next set of vertices
                for next_vertex in next_vertices_relation(graph, curr_node) #for dominator this is the children of the current node, unconstrained. for postdominator it is the parents, unconstrained.
                    if !(next_vertex in visited)
                        push!(work_heap, next_vertex)
                        push!(visited, next_vertex)
                    end

                end
            end
        end
    end

    return current_dom
end

"""Computes immediate dominators/postdominators of all nodes in a graph. Assumes `predecessors` contains the predecessors of nodes in post order for dominators (the leaves of the graph have the lowest post order numbers and come earliest in predecessors), and reverse post order for postdominators (similar except order of nodes reverse). In the immediate **dominator** case the root node is visited first and is its own dominator. The remaining nodes are visited in reverse post order. The immediate dominator of any node i = node i ∪ (∩ dominance(predecessors[i]))."""
function simple_dominance(predecessors::Vector{Vector{Int64}}, dominance::Union{Nothing,Vector{BitVector}}=nothing, idoms::Union{Nothing,Vector{Int64}}=nothing)
    if dominance === nothing
        dominance = [falses(length(predecessors)) for _ in 1:length(predecessors)]
    end

    temp = BitVector(undef, length(predecessors))
    if idoms === nothing
        idoms = Vector{Int64}(undef, length(predecessors))
    end

    for index in length(dominance):-1:1
        dominance[index][index] = 1
        if index != length(dominance)
            temp .= dominance[predecessors[index][1]] #first value in intersection of dom_masks of predecessors. This intersection contains all the dominators of this node
        else
            temp .= dominance[index] #root or varible node gets its own dom_mask
        end

        for i in 2:min(2, length(predecessors[index]))
            @. temp = temp & dominance[predecessors[index][i]]
        end
        @. dominance[index] |= temp
    end

    for (i, mask) in pairs(dominance)
        if i ≠ length(predecessors) #don't reset root or variable node mask since it is the only node that can be its own idom.
            mask[i] = 0 #set dom mask for this node index to 0 so don't incorrectly find it to be its own idom.
            idoms[i] = findfirst(mask) #nodes are in reverse post order for doms and post order for pdoms. The index of the first non-zero bit in mask is the immediate dominator of node i.
        else
            idoms[i] = i
        end
    end

    return idoms
end
export simple_dominance


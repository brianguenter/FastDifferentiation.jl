abstract type AbstractFactorableSubgraph end
export AbstractFactorableSubgraph
abstract type DominatorSubgraph <: AbstractFactorableSubgraph end
abstract type PostDominatorSubgraph <: AbstractFactorableSubgraph end

struct FactorableSubgraph{T<:Integer,S<:AbstractFactorableSubgraph}
    graph::DerivativeGraph{T}
    subgraph::Tuple{T,T}
    times_used::T
    reachable_roots::BitVector
    reachable_variables::BitVector
    dom_mask::Union{Nothing,BitVector}
    pdom_mask::Union{Nothing,BitVector}
    #if these two numbers are unchanged since the creation of the subgraph then the subgraph still exists. Quick test for subgraph destruction.
    num_dominator_edges::T #number of edges from the dominator node that satisfy the subgraph relation.
    num_dominated_edges::T #number of edges from the dominated node that satisfy the subgraph relation

    function FactorableSubgraph{T,DominatorSubgraph}(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, dom_mask::BitVector, roots_reachable::BitVector, variables_reachable::BitVector) where {T<:Integer}
        @assert dominating_node > dominated_node
        dominator_relation_edges = length(child_edges(graph, dominating_node)) #when graph is created all children of dominating_node satisfy relation_edges. After subgraphs are factored this may no longer be true. 
        constraint = PathConstraint(dominating_node, graph, true, dom_mask, variables_reachable)
        tmp_edges = get_edge_vector()
        relation_edges!(constraint, dominated_node, tmp_edges)
        dominated_relation_edges = length(tmp_edges)
        reclaim_edge_vector(tmp_edges)
        sum(dom_mask) * sum(variables_reachable)
        return new{T,DominatorSubgraph}(graph, (dominating_node, dominated_node), sum(dom_mask) * sum(variables_reachable), roots_reachable, variables_reachable, dom_mask, nothing, dominator_relation_edges, dominated_relation_edges)
    end

    function FactorableSubgraph{T,PostDominatorSubgraph}(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, pdom_mask::BitVector, roots_reachable::BitVector, variables_reachable::BitVector) where {T<:Integer}
        @assert dominating_node < dominated_node
        dominator_relation_edges = length(parent_edges(graph, dominating_node)) #when graph is created all parents of dominating_node satisfy relation_edges. After subgraphs are factored this may no longer be true. 
        constraint = PathConstraint(dominating_node, graph, false, roots_reachable, pdom_mask)
        tmp_edges = get_edge_vector()
        relation_edges!(constraint, dominated_node, tmp_edges)
        dominated_relation_edges = length(tmp_edges)
        reclaim_edge_vector(tmp_edges)
        return new{T,PostDominatorSubgraph}(graph, (dominating_node, dominated_node), sum(roots_reachable) * sum(pdom_mask), roots_reachable, variables_reachable, nothing, pdom_mask, dominator_relation_edges, dominated_relation_edges)
    end
end
export FactorableSubgraph

FactorableSubgraph(args::Tuple) = FactorableSubgraph(args...)

dominator_subgraph(args) = dominator_subgraph(args...)
dominator_subgraph(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, dom_mask::BitVector, roots_reachable::BitVector, variables_reachable::BitVector) where {T<:Integer} = FactorableSubgraph{T,DominatorSubgraph}(graph, dominating_node, dominated_node, dom_mask, roots_reachable, variables_reachable)
dominator_subgraph(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, dom_mask::S, roots_reachable::S, variables_reachable::S) where {T<:Integer,S<:Vector{Bool}} = dominator_subgraph(graph, dominating_node, dominated_node, BitVector(dom_mask), BitVector(roots_reachable), BitVector(variables_reachable))
export dominator_subgraph

postdominator_subgraph(args) = postdominator_subgraph(args...)
postdominator_subgraph(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, pdom_mask::BitVector, roots_reachable::BitVector, variables_reachable::BitVector) where {T<:Integer} = FactorableSubgraph{T,PostDominatorSubgraph}(graph, dominating_node, dominated_node, pdom_mask, roots_reachable, variables_reachable)
postdominator_subgraph(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, pdom_mask::S, roots_reachable::S, variables_reachable::S) where {T<:Integer,S<:Vector{Bool}} = postdominator_subgraph(graph, dominating_node, dominated_node, BitVector(pdom_mask), BitVector(roots_reachable), BitVector(variables_reachable))
export postdominator_subgraph

graph(a::FactorableSubgraph) = a.graph

"""Returns a tuple of ints (dominator vertex,dominated vertex) that are the top and bottom vertices of the subgraph"""
vertices(subgraph::FactorableSubgraph) = subgraph.subgraph
export vertices

reachable_variables(a::FactorableSubgraph) = a.reachable_variables
function mask_variables!(a::FactorableSubgraph, mask::BitVector)
    @assert domain_dimension(graph(a)) == length(mask)
    a.reachable_variables .&= mask
end

reachable_roots(a::FactorableSubgraph) = a.reachable_roots
function mask_roots!(a::FactorableSubgraph, mask::BitVector)
    @assert codomain_dimension(graph(a)) == length(mask)
    a.reachable_roots .&= mask
end

reachable(a::FactorableSubgraph{T,DominatorSubgraph}) where {T} = reachable_variables(a)
reachable(a::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = reachable_roots(a)
export reachable

# """Returns reachable roots of the node at `node_index`"""
# function reachable(a::FactorableSubgraph{T,DominatorSubgraph}, node_index::T) where {T}
#     gr = graph(a)
#     result = falses(codomain_dimension(a))

#     for edge in parent_edges(graph, node_index)
#         @. result |= reachable_roots(edge)
#     end
#     return result
# end

# """Returns reachable variables of the node at `node_index`"""
# function reachable(a::FactorableSubgraph{T,DominatorSubgraph}, node_index::T) where {T}
#     gr = graph(a)
#     result = falses(codomain_dimension(a))

#     for edge in child_edges(graph, node_index)
#         @. result |= reachable_variables(edge)
#     end
#     return result
# end

reachable_dominance(a::FactorableSubgraph{T,DominatorSubgraph}) where {T} = a.dom_mask
reachable_dominance(a::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = a.pdom_mask
export reachable_dominance

dominating_node(a::FactorableSubgraph{T,S}) where {T,S<:Union{DominatorSubgraph,PostDominatorSubgraph}} = a.subgraph[1]
export dominating_node
dominated_node(a::FactorableSubgraph{T,S}) where {T,S<:Union{DominatorSubgraph,PostDominatorSubgraph}} = a.subgraph[2]
export dominated_node

times_used(a::FactorableSubgraph) = a.times_used
export times_used

node_difference(a::FactorableSubgraph) = abs(a.subgraph[1] - a.subgraph[2])

function Base.show(io::IO, a::FactorableSubgraph)
    print(io, summarize(a))
end

function summarize(a::FactorableSubgraph{T,DominatorSubgraph}) where {T}
    doms = ""
    doms *= to_string(reachable_dominance(a), "r")
    doms *= " ↔ "
    doms *= to_string(reachable_variables(a), "v")

    return "[" * doms * " $(times_used(a))* " * string((vertices(a))) * "]"
end
export summarize

function summarize(a::FactorableSubgraph{T,PostDominatorSubgraph}) where {T}
    doms = ""
    doms *= to_string(reachable_roots(a), "r")
    doms *= " ↔ "
    doms *= to_string(reachable_dominance(a), "v")

    return "[" * doms * " $(times_used(a))* " * string((vertices(a))) * "]"
end
export summarize

"""Returns parent edges if subgraph is dominator and child edges otherwise. Parent edges correspond to the forward traversal of a dominator subgraph in graph factorization, analogously for postdominator subgraph"""
forward_edges(a::FactorableSubgraph{T,DominatorSubgraph}, edge::PathEdge) where {T} = parent_edges(graph(a), edge)
forward_edges(a::FactorableSubgraph{T,PostDominatorSubgraph}, edge::PathEdge) where {T} = child_edges(graph(a), edge)

forward_edges(a::FactorableSubgraph{T,DominatorSubgraph}, node_index::T) where {T} = parent_edges(graph(a), node_index)
forward_edges(a::FactorableSubgraph{T,PostDominatorSubgraph}, node_index::T) where {T} = child_edges(graph(a), node_index)
export forward_edges

"""Returns child edges if subgraph is dominator and parent edges otherwise. Child edges correspond to the backward check for paths bypassing the dominated node of a dominator subgraph, analogously for postdominator subgraph"""
backward_edges(a::FactorableSubgraph{T,DominatorSubgraph}, node_index::T) where {T} = child_edges(graph(a), node_index)
backward_edges(a::FactorableSubgraph{T,PostDominatorSubgraph}, node_index::T) where {T} = parent_edges(graph(a), node_index)

backward_edges(a::FactorableSubgraph, edge::PathEdge) = backward_edges(a, forward_vertex(a, edge))

test_edge(a::FactorableSubgraph{T,DominatorSubgraph}, edge::PathEdge) where {T} = subset(reachable_dominance(a), reachable_roots(edge)) && overlap(reachable_variables(a), reachable_variables(edge))
test_edge(a::FactorableSubgraph{T,PostDominatorSubgraph}, edge::PathEdge) where {T} = subset(reachable_dominance(a), reachable_variables(edge)) && overlap(reachable_roots(a), reachable_roots(edge))

reachable_dominance(::FactorableSubgraph{T,DominatorSubgraph}, edge::PathEdge) where {T} = reachable_roots(edge)
reachable_dominance(::FactorableSubgraph{T,PostDominatorSubgraph}, edge::PathEdge) where {T} = reachable_variables(edge)

non_dominance_mask(::FactorableSubgraph{T,DominatorSubgraph}, edge::PathEdge) where {T} = reachable_variables(edge)
non_dominance_mask(::FactorableSubgraph{T,PostDominatorSubgraph}, edge::PathEdge) where {T} = reachable_roots(edge)

non_dominance_mask(a::FactorableSubgraph{T,DominatorSubgraph}) where {T} = reachable_variables(a)
non_dominance_mask(a::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = reachable_roots(a)

non_dominance_dimension(subgraph::FactorableSubgraph{T,DominatorSubgraph}) where {T} = domain_dimension(graph(subgraph))
non_dominance_dimension(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = codomain_dimension(graph(subgraph))

forward_vertex(::FactorableSubgraph{T,DominatorSubgraph}, edge::PathEdge) where {T} = top_vertex(edge)
forward_vertex(::FactorableSubgraph{T,PostDominatorSubgraph}, edge::PathEdge) where {T} = bott_vertex(edge)

function next_valid_edge(a::FactorableSubgraph, current_edge::PathEdge{T}) where {T}
    if forward_vertex(a, current_edge) == dominating_node(a) #reached the end of the subgraph
        return nothing
    else
        local edge_next::PathEdge{T}
        count = 0
        for edge in forward_edges(a, current_edge) #should always be a next edge because top_vertex(current_edge) != dominance_node(a)
            if test_edge(a, edge)
                count += 1
                # @assert count ≤ 1 #in a properly processed subgraph there should not be branches on paths from dominated to dominating node.
                if count > 1
                    return nothing
                end
                edge_next = edge
            end
        end
        if count == 0
            return nothing #no valid next edge. This can occur if reachable variables or roots were reset by a previous factorization step.
        else
            return edge_next
        end
    end
end

function isa_connected_path(a::FactorableSubgraph, start_edge::PathEdge{T}) where {T}
    current_edge = start_edge

    if test_edge(a, start_edge) #ensure that start_edge satisfies conditions for being on a connected path
        while forward_vertex(a, current_edge) != dominating_node(a)
            current_edge = next_valid_edge(a, current_edge)
            if current_edge === nothing
                return false
            end
        end
        return true
    else
        return false
    end
end
export isa_connected_path

function edges_on_path(a::FactorableSubgraph, start_edge::PathEdge{T}) where {T}
    current_edge = start_edge
    result = PathEdge{T}[]

    if test_edge(a, start_edge) #ensure that start_edge satisfies conditions for being on a connected path
        while forward_vertex(a, current_edge) != dominating_node(a)
            current_edge = next_valid_edge(a, current_edge)
            if current_edge === nothing
                return (false, PathEdge{T}[])
            else
                push!(result, current_edge)
            end
        end
        return (true, result)
    else
        return (false, PathEdge{T}[])
    end
end
export edges_on_path

"""Splits edges which have roots not in the `dominance_mask` of `subgraph`. Original edge has only roots in `dominance_mask`. A new edge is added to the graph that contains only roots not in `dominance_mask`."""
function add_non_dom_edges!(subgraph::FactorableSubgraph{T,S}) where {T,S<:AbstractFactorableSubgraph}
    temp_edges = PathEdge{T}[]

    for s_edge in forward_edges(subgraph, dominated_node(subgraph))
        if test_edge(subgraph, s_edge)
            for pedge in edge_path(subgraph, s_edge)
                edge_mask = reachable_dominance(subgraph, pedge)
                diff = set_diff(edge_mask, reachable_dominance(subgraph)) #important that diff is a new BitVector, not reused.
                if any(diff)
                    gr = graph(subgraph)

                    if S === DominatorSubgraph
                        push!(temp_edges, PathEdge(top_vertex(pedge), bott_vertex(pedge), value(pedge), copy(reachable_variables(pedge)), diff)) #create a new edge that accounts for roots not in the dominance mask
                    else

                        push!(temp_edges, PathEdge(top_vertex(pedge), bott_vertex(pedge), value(pedge), diff, copy(reachable_roots(pedge)))) #create a new edge that accounts for roots not in the     dominance mask    
                    end

                    @. edge_mask &= !diff #in the original edge reset the roots/variables not in dominance mask
                end
            end
        end
    end
    gr = graph(subgraph)
    for edge in temp_edges
        add_edge!(gr, edge)
    end

    return nothing
end
export add_non_dom_edges!

"""Sets the reachable root and variable masks for every edge in `DominatorSubgraph` `subgraph`. """
function reset_edge_masks!(subgraph::FactorableSubgraph{T}) where {T}
    edges_to_delete = PathEdge{T}[]

    for fwd_edge in forward_edges(subgraph, dominated_node(subgraph))
        bypass_mask = .!copy(non_dominance_mask(subgraph)) #bypass mask tracks which variables/roots are on a backward path that bypasses the dominated_node. These variables/roots cannot be reset. 0 means can be reset 1 means can't.

        if test_edge(subgraph, fwd_edge)
            for pedge in edge_path(subgraph, fwd_edge)

                mask = non_dominance_mask(subgraph, pedge)
                @. mask = mask & bypass_mask #if any bits in bypass mask are 1 then those bits won't be reset. 

                if !any(bypass_mask .& non_dominance_mask(subgraph)) #no edges bypass dominated node so can reset all dominance bits. If any edges bypass cannot reset any dominance bits.
                    fmask = reachable_dominance(subgraph, fwd_edge)
                    fmask = fmask .& .!reachable_dominance(subgraph)
                end

                if forward_vertex(subgraph, pedge) != dominating_node(subgraph)
                    for bedge in backward_edges(subgraph, pedge)
                        if test_edge(subgraph, bedge) && bedge !== pedge
                            #want to test by object identity - don't want to include the non_dom mask of the current edge 
                            bypass_mask .|= non_dominance_mask(subgraph, bedge)
                        end
                    end
                end

                if can_delete(pedge)
                    push!(edges_to_delete, pedge)
                end
            end
        end
    end

    return edges_to_delete
end
export reset_edge_masks!


function check_edges(subgraph::FactorableSubgraph, edge_list::Vector{PathEdge{T}}) where {T}
    #make sure have at least two edges that are on a valid path from dominated to dominating node
    count = 0
    for edge in edge_list
        if test_edge(subgraph, edge)
            count += 1
        end
    end
    if count < 2
        return false
    else
        return true
    end
end

"""Returns all edges in the subgraph as a set or `nothing` if the subgraph no longer exists."""
function subgraph_edges(subgraph::FactorableSubgraph{T}) where {T}
    result = Set{PathEdge{T}}()

    if !subgraph_exists(subgraph)
        return nothing
    else
        for fedge in forward_edges(subgraph, dominated_node(subgraph))
            good_edges, tmp = edges_on_path(subgraph, fedge)
            #this needs to be recursive. Can't use edges_on_path because that assumes the subgraph doesn't have branches. subgraph_edges is only called when branches have edges.
            if good_edges
                union!(result, tmp)
                push!(result, fedge)
            end
        end

        @assert length(result) ≥ 2 "If subgraph exists should be at least two valid edges in subgraph. Instead got none."

        return result
    end
end
export subgraph_edges


"""returns subgraph edges, as a `Set`, and nodes, as a `Vector`."""
function deconstruct_subgraph(subgraph::FactorableSubgraph)
    sub_edges = subgraph_edges(subgraph)
    sub_nodes = map(x -> bott_vertex(x), collect(sub_edges))
    if subgraph isa FactorableSubgraph{T,DominatorSubgraph} where {T}
        push!(sub_nodes, dominating_node(subgraph))
    else
        push!(sub_nodes, dominated_node(subgraph))
    end

    return sub_edges, unique(sub_nodes)
end
export deconstruct_subgraph

"""Returns true if the subgraph is still a factorable dominance subgraph, false otherwise"""
function subgraph_exists(subgraph::FactorableSubgraph)
    #Do fast tests that guarantee subgraph has been destroyed by factorization: no edges connected to dominated node, dominated_node or dominator node has < 2 subgraph edges
    #This is inefficient since many tests require just the number of edges but this code creates temp arrays containing the edges and then measures the length. Optimize later by having separate children and parents fields in edges structure of RnToRmGraph. Then num_parents and num_children become fast and allocation free.

    #need at least two parent edges from dominated_node or subgraph doesn't exist
    fedges = forward_edges(subgraph, dominated_node(subgraph))
    bedges = backward_edges(subgraph, dominating_node(subgraph))

    if length(fedges) < 2 || length(bedges) < 2 #need at least two forward edges from dominated_node and two backward edges from dominating node or subgraph doesn't exist
        return false
    elseif !check_edges(subgraph, fedges) #verify that all edges have correct reachability to be on a valid path from dominated node to dominating node
        return false
    elseif !check_edges(subgraph, bedges)
        return false
    else
        count = 0
        sub_edges = Set{PathEdge}()

        for edge in fedges
            if isa_connected_path(subgraph, edge) !== nothing
                count += 1
            end

            @invariant begin
                good_edges, tmp = edges_on_path(subgraph, edge)
                good_subgraph = true

                if good_edges
                    for pedge in tmp
                        if in(pedge, sub_edges)
                            good_subgraph = false
                            break
                        end
                        push!(sub_edges, pedge)
                    end
                end
                good_subgraph
            end "Visited edge $(vertices(pedge)) more than once in subgraph $(vertices(subgraph)). Edges in factorable subgraph should only be visited once. "
        end
        if count >= 2
            return true
        else
            return false
        end
    end
end

#PathIterator has redundant computation because path is checked for connectivity when creating iterator and then path is traversed again when running iterator. Not sure how if it is possible to use Iterator framework without doing this, or making the calling code much more complex.
struct PathIterator{T<:Integer,S<:FactorableSubgraph}
    subgraph::S
    start_edge::PathEdge{T}

    function PathIterator(subgraph::S, start_edge::PathEdge{T}) where {T,S<:FactorableSubgraph{T,DominatorSubgraph}}
        @assert bott_vertex(start_edge) == dominated_node(subgraph)
        return new{T,S}(subgraph, start_edge)
    end

    function PathIterator(subgraph::S, start_edge::PathEdge{T}) where {T,S<:FactorableSubgraph{T,PostDominatorSubgraph}}
        @assert top_vertex(start_edge) == dominated_node(subgraph)
        return new{T,S}(subgraph, start_edge)
    end
end
export PathIterator

edge_path(subgraph::FactorableSubgraph, start_edge) = PathIterator(subgraph, start_edge)
export edge_path

"""Returns an iterator for a single path in a factorable subgraph. If the path has been destroyed by factorization returns nothing."""
function Base.iterate(a::PathIterator{T,S}) where {T,S<:FactorableSubgraph}
    if !isa_connected_path(a.subgraph, a.start_edge)
        return nothing
    else
        return (a.start_edge, a.start_edge)
    end
end

function Base.iterate(a::PathIterator{T,S}, state::PathEdge{T}) where {T,S<:FactorableSubgraph}
    if forward_vertex(a.subgraph, state) == dominating_node(a.subgraph)
        return nothing
    else
        edge_next = next_valid_edge(a.subgraph, state)

        @assert edge_next !== nothing #tested for connected path when creating iterator so should never get nothing return because edge is not at the dominator node

        return (edge_next, edge_next)
    end
end

Base.IteratorSize(::Type{<:PathIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:PathIterator}) = Base.HasEltype()
Base.eltype(::Type{<:PathIterator{T}}) where {T} = PathEdge{T}


function r1r1subgraph_edges(graph::DerivativeGraph{T}, root_index::T, variable_index::T) where {T}
    visited = Set{T}()
    sub_edges = Set{PathEdge}()
    return _r1r1subgraph_edges(graph, root_index_to_postorder_number(graph, root_index), root_index, variable_index, visited, sub_edges)
end


"""Returns the edges in the subgraph connecting `root_index` and `variable_index`. This is an R¹->R¹ function. Used for debugging."""
function _r1r1subgraph_edges(graph::DerivativeGraph{T},
    current_index::T,
    root_index::T,
    variable_index::T,
    visited::Set{T},
    sub_edges::Set{PathEdge}) where {T}

    if !in(current_index, visited)
        push!(visited, current_index)
        for child in child_edges(graph, current_index)
            if reachable_roots(child)[root_index]#&& reachable_variables(child)[variable_index]
                push!(sub_edges, child)
                _r1r1subgraph_edges(graph, bott_vertex(child), root_index, variable_index, visited, sub_edges)
            end
        end
    end

    return sub_edges
end

"""Finds the node indices in the subgraph. Computes relations for each node (indices of parents inside subgraph for dom, children for pdom. Returns two vectors: the indices of the graph nodes, and a vector of vectors that contains the relations informations."""
function graph_array(subgraph::DominatorSubgraph)

end




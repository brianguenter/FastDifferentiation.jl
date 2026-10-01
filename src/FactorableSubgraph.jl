abstract type AbstractFactorableSubgraph end
abstract type DominatorSubgraph <: AbstractFactorableSubgraph end
abstract type PostDominatorSubgraph <: AbstractFactorableSubgraph end

"""
    FactorableSubgraph{T<:Integer,S<:AbstractFactorableSubgraph}

# Fields

- `graph::DerivativeGraph{T}`
- `subgraph::Tuple{T,T}`
- `times_used::T`
- `reachable_roots::BitVector`
- `reachable_variables::BitVector`
- `dom_mask::Union{Nothing,BitVector}`
- `pdom_mask::Union{Nothing,BitVector}`
"""
struct FactorableSubgraph{T<:Integer,S<:AbstractFactorableSubgraph}
    graph::DerivativeGraph{T}
    subgraph::Tuple{T,T}
    times_used::T
    reachable_roots::BitVector
    reachable_variables::BitVector
    dom_mask::Union{Nothing,BitVector}
    pdom_mask::Union{Nothing,BitVector}

    function FactorableSubgraph{T,DominatorSubgraph}(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, dom_mask::BitVector, roots_reachable::BitVector, variables_reachable::BitVector) where {T<:Integer}
        @assert dominating_node > dominated_node

        return new{T,DominatorSubgraph}(graph, (dominating_node, dominated_node), sum(dom_mask) * sum(variables_reachable), roots_reachable, variables_reachable, dom_mask, nothing)
    end

    function FactorableSubgraph{T,PostDominatorSubgraph}(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, pdom_mask::BitVector, roots_reachable::BitVector, variables_reachable::BitVector) where {T<:Integer}
        @assert dominating_node < dominated_node

        return new{T,PostDominatorSubgraph}(graph, (dominating_node, dominated_node), sum(roots_reachable) * sum(pdom_mask), roots_reachable, variables_reachable, nothing, pdom_mask)
    end
end

dominator_subgraph(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, dom_mask::BitVector, roots_reachable::BitVector, variables_reachable::BitVector) where {T<:Integer} = FactorableSubgraph{T,DominatorSubgraph}(graph, dominating_node, dominated_node, dom_mask, roots_reachable, variables_reachable)
dominator_subgraph(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, dom_mask::S, roots_reachable::S, variables_reachable::S) where {T<:Integer,S<:Vector{Bool}} = dominator_subgraph(graph, dominating_node, dominated_node, BitVector(dom_mask), BitVector(roots_reachable), BitVector(variables_reachable))


postdominator_subgraph(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, pdom_mask::BitVector, roots_reachable::BitVector, variables_reachable::BitVector) where {T<:Integer} = FactorableSubgraph{T,PostDominatorSubgraph}(graph, dominating_node, dominated_node, pdom_mask, roots_reachable, variables_reachable)
postdominator_subgraph(graph::DerivativeGraph{T}, dominating_node::T, dominated_node::T, pdom_mask::S, roots_reachable::S, variables_reachable::S) where {T<:Integer,S<:Vector{Bool}} = postdominator_subgraph(graph, dominating_node, dominated_node, BitVector(pdom_mask), BitVector(roots_reachable), BitVector(variables_reachable))


graph(a::FactorableSubgraph) = a.graph

"""
    vertices(subgraph::FactorableSubgraph)

Returns a tuple of ints (dominator vertex,dominated vertex) that are the top and bottom vertices of the subgraph"""
vertices(subgraph::FactorableSubgraph) = subgraph.subgraph


reachable_variables(a::FactorableSubgraph) = a.reachable_variables

reachable_roots(a::FactorableSubgraph) = a.reachable_roots

reachable_dominance(a::FactorableSubgraph{T,DominatorSubgraph}) where {T} = a.dom_mask
reachable_dominance(a::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = a.pdom_mask


dominating_node(a::FactorableSubgraph{T,S}) where {T,S<:Union{DominatorSubgraph,PostDominatorSubgraph}} = a.subgraph[1]

dominated_node(a::FactorableSubgraph{T,S}) where {T,S<:Union{DominatorSubgraph,PostDominatorSubgraph}} = a.subgraph[2]


times_used(a::FactorableSubgraph) = a.times_used


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


function summarize(a::FactorableSubgraph{T,PostDominatorSubgraph}) where {T}
    doms = ""
    doms *= to_string(reachable_roots(a), "r")
    doms *= " ↔ "
    doms *= to_string(reachable_dominance(a), "v")

    return "[" * doms * " $(times_used(a))* " * string((vertices(a))) * "]"
end


"""
    forward_edges(a::FactorableSubgraph{T,DominatorSubgraph}, node_index::T)

Returns parent edges if subgraph is dominator and child edges otherwise. Parent edges correspond to the forward traversal of a dominator subgraph in graph factorization, analogously for postdominator subgraph"""
forward_edges(a::FactorableSubgraph{T,DominatorSubgraph}, node_index::T) where {T} = parent_edges(graph(a), node_index)
forward_edges(a::FactorableSubgraph{T,PostDominatorSubgraph}, node_index::T) where {T} = child_edges(graph(a), node_index)


"""
    backward_edges(a::FactorableSubgraph{T,DominatorSubgraph}, node_index::T)

Returns child edges if subgraph is dominator and parent edges otherwise. Child edges correspond to the backward check for paths bypassing the dominated node of a dominator subgraph, analogously for postdominator subgraph"""
backward_edges(a::FactorableSubgraph{T,DominatorSubgraph}, node_index::T) where {T} = child_edges(graph(a), node_index)
backward_edges(a::FactorableSubgraph{T,PostDominatorSubgraph}, node_index::T) where {T} = parent_edges(graph(a), node_index)

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

"""
    add_non_dom_edges!(subgraph::FactorableSubgraph{T,S})

Splits boundary edges leaving `dominated_node(subgraph)` that have roots/variables not in the `dominance_mask` of `subgraph`. The original edge retains only roots/variables in `dominance_mask`. A new edge is added to the graph that contains only roots/variables not in `dominance_mask`."""
function add_non_dom_edges!(subgraph::FactorableSubgraph{T,S}) where {T,S<:AbstractFactorableSubgraph}
    temp_edges = PathEdge{T}[]

    for s_edge in forward_edges(subgraph, dominated_node(subgraph))
        if test_edge(subgraph, s_edge)
            edge_mask = reachable_dominance(subgraph, s_edge)
            diff = set_diff(edge_mask, reachable_dominance(subgraph)) #important that diff is a new BitVector, not reused.
            if any(diff)
                if S === DominatorSubgraph
                    push!(temp_edges, PathEdge(top_vertex(s_edge), bott_vertex(s_edge), value(s_edge), copy(reachable_variables(s_edge)), diff)) #create a new edge that accounts for roots not in the dominance mask
                else
                    push!(temp_edges, PathEdge(top_vertex(s_edge), bott_vertex(s_edge), value(s_edge), diff, copy(reachable_roots(s_edge)))) #create a new edge that accounts for roots not in the dominance mask
                end

                @. edge_mask &= !diff #in the original edge reset the roots/variables not in dominance mask
            end
        end
    end
    gr = graph(subgraph)
    for edge in temp_edges
        add_edge!(gr, edge)
    end

    return nothing
end


"""
    reset_edge_masks!(subgraph::FactorableSubgraph{T})

Resets the non-dominance reachable masks for the boundary edges leaving `dominated_node(subgraph)` that are factored by this partition, returning any edges that can now be deleted."""
function reset_edge_masks!(subgraph::FactorableSubgraph{T}) where {T}
    edges_to_delete = PathEdge{T}[]
    bypass_mask = .!copy(non_dominance_mask(subgraph))

    for fwd_edge in forward_edges(subgraph, dominated_node(subgraph))
        if test_edge(subgraph, fwd_edge)
            mask = non_dominance_mask(subgraph, fwd_edge)
            @. mask = mask & bypass_mask

            if can_delete(fwd_edge)
                push!(edges_to_delete, fwd_edge)
            end
        end
    end

    return edges_to_delete
end

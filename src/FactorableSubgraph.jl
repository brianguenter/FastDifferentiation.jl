abstract type AbstractFactorableSubgraph end
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

subgraph_vertices(subgraph::FactorableSubgraph) = subgraph.subgraph
export subgraph_vertices

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

dominance_mask(a::FactorableSubgraph{T,DominatorSubgraph}) where {T} = a.dom_mask
dominance_mask(a::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = a.pdom_mask
export dominance_mask

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
    doms *= to_string(dominance_mask(a), "r")
    doms *= " ↔ "
    doms *= to_string(reachable_variables(a), "v")

    return "[" * doms * " $(times_used(a))* " * string((subgraph_vertices(a))) * "]"
end
export summarize

function summarize(a::FactorableSubgraph{T,PostDominatorSubgraph}) where {T}
    doms = ""
    doms *= to_string(reachable_roots(a), "r")
    doms *= " ↔ "
    doms *= to_string(dominance_mask(a), "v")

    return "[" * doms * " $(times_used(a))* " * string((subgraph_vertices(a))) * "]"
end
export summarize

#TODO delete this if still not used
function subgraph_nodes(subgraph::FactorableSubgraph{T,DominatorSubgraph}) where {T}
    result = T[]

    start_node = dominated_node(subgraph)
    constraint = x -> subset(dominance_mask(subgraph), reachable_roots(x)) && overlap(reachable_variables(subgraph), reachable_variables(x))
    push!(result, start_node)
    current_level = Set{T}(start_node)
    next_level = Set{T}()
    new_nodes = true

    while new_nodes
        new_nodes = false
        for begin_node in current_level
            for next_node in ConstrainedPathIterator(begin_node, dominating_node(subgraph), node_edges(graph(subgraph), begin_node), true, constraint)
                if !in(next_node, next_level)
                    push!(next_level, next_node)
                    new_nodes = true
                end
            end
        end
        append!(result, values(next_level))
    end

    return sort!(result) #return sorted by postorder number with largest postorder number at the end of the array
end


"""Traverse edges in `subgraph` to see if `dominating_node(subgraph)` is still the idom of `dominated_node(subgraph)` or if factorization has destroyed the subgraph. If there are two or more paths from dominated to dominating node and there are no branches on these paths then the subgraph still exists."""
function valid_paths(constraint, subgraph::FactorableSubgraph{T,S}) where {T,S<:Union{PostDominatorSubgraph,DominatorSubgraph}}
    start_edges = get_edge_vector()
    relation_edges!(constraint, dominated_node(subgraph), start_edges)

    if length(start_edges) == 0 || length(start_edges) == 1
        reclaim_edge_vector(start_edges)
        return false #subgraph has been destroyed
    else
        count = 0
        for pedge in start_edges
            path_edges = get_edge_vector()
            flag = edges_on_path!(constraint, dominating_node(subgraph), S == DominatorSubgraph, pedge, path_edges)
            if flag == 1
                count += 1
            elseif flag == 2 #there is a branch in the edge path which can only happen if the subgraph has been destroyed
                reclaim_edge_vector(path_edges)
                reclaim_edge_vector(start_edges)
                return false
            end
        end

        if count < 2
            reclaim_edge_vector(start_edges)
            return false #don't have two non-branching paths from dominated to dominating node so subgraph has been destroyed
        end
    end
    reclaim_edge_vector(start_edges)
    return true
end

"""Sets the reachable root and variable masks for every edge in `DominatorSubgraph` `subgraph`. """
function reset_edge_masks!(subgraph::FactorableSubgraph{T,DominatorSubgraph}) where {T}
    gr = graph(subgraph)
    dominator = dominating_node(subgraph)
    dominated = dominated_node(subgraph)
    edge_constraint = next_edge_constraint(subgraph)


    for start_edge in relation_edges!(edge_constraint, dominated)
        current_edge = start_edge
        bypass_mask = trues(codomain_dimension(gr)) .& .!reachable_variables(subgraph) #if bypassmask is 0 at index i then variable vᵢ can be reset. Initially all reachable variables of the subgraph are resettable.

        while true
            current_node = top_vertex(current_edge)
            rdiff = set_diff(reachable_roots(current_edge), reachable_roots(subgraph))
            if is_zero(rdiff)
                rvars = reachable_variables(current_edge)
                rvars .= reachable_variables(current_edge) .& bypass_mask
            end

            if is_zero(bypass_mask .& reachable_variables(subgraph)) #no child edge paths bypass dominated node so can reset all rᵢ ∈ Rdom.
                rroots = reachable_roots(current_edge)
                rroots .= reachable_roots(current_edge) .& .!dominance_mask(subgraph)
            end

            if is_zero(reachable_variables(current_edge)) || is_zero(reachable_roots(current_edge)) #paths to roots or variables so cannot participate in any future path products. Remove edge to simplify graph.
                delete_edge!(gr, current_edge)
            end

            if current_node == dominator
                break
            end

            for child in child_edges(gr, current_node)
                if child !== current_edge
                    bypass_mask .|= reachable_variables(child)
                end
            end
        end
    end
end


"""Sets the reachable root and variable masks for every edge in `PostDominatorSubgraph` `subgraph`. """
function reset_edge_masks!(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}) where {T}
    gr = graph(subgraph)
    dominator = dominating_node(subgraph)
    dominated = dominated_node(subgraph)
    edge_constraint = next_edge_constraint(subgraph)


    for start_edge in relation_edges!(edge_constraint, dominated)
        current_edge = start_edge
        bypass_mask = trues(domain_dimension(gr)) .& .!reachable_roots(subgraph) #if bypassmask is 0 at index i then variable vᵢ can be reset. Initially all reachable variables of the subgraph are resettable.

        while true
            current_node = bott_vertex(current_edge)
            vdiff = set_diff(reachable_variables(current_edge), reachable_variables(subgraph))
            if is_zero(vdiff)
                rroots = reachable_roots(current_edge)
                rrots .= reachable_roots(current_edge) .& bypass_mask
            end

            if is_zero(bypass_mask .& reachable_roots(subgraph)) #no child edge paths bypass dominated node so can reset all rᵢ ∈ Rdom.
                rvars = reachable_variables(current_edge)
                rvars .= reachable_variables(current_edge) .& .!dominance_mask(subgraph)
            end

            if current_node == dominator
                break
            end

            for child in parent_edges(gr, current_node)
                if child !== current_edge
                    bypass_mask .|= reachable_roots(child)
                end
            end
        end
    end
end

"""Returns true if the subgraph is still a factorable dominance subgraph, false otherwise"""
function subgraph_exists(subgraph::FactorableSubgraph{T,DominatorSubgraph}) where {T}
    #Do fast tests that guarantee subgraph has been destroyed by factorization: no edges connected to dominated node, dominated_node or dominator node has < 2 subgraph edges
    #This is inefficient since many tests require just the number of edges but this code creates temp arrays containing the edges and then measures the length. Optimize later by having separate children and parents fields in edges structure of RnToRmGraph. Then num_parents and num_children become fast and allocation free.

    constraint = next_edge_constraint(subgraph)
    dgraph = graph(subgraph)

    sub_edges = get_edge_vector()
    relation_edges!(constraint, dominated_node(subgraph), sub_edges) #need at least two parent edges from dominated_node or subgraph doesn't exist
    if sub_edges === nothing
        return false
    elseif length(sub_edges) < 2
        return false
    else
        cedges = child_edges(dgraph, dominating_node(subgraph)) #need at least two unconstrained child edges from dominating_node or subgraph doesn't exist
        if length(cedges) < 2
            return false
        else
            edge_count = 0

            for x in cedges

                #not certain this test is correct. Only correct if subset(dominance_mask(subgraph), reachable_roots(x)) == true means there are no root values for which the subgraph is factorable. But this will only be true if overlap(dominance_mask)
                if subset(dominance_mask(subgraph), reachable_roots(x)) && overlap(reachable_variables(subgraph), reachable_variables(x))
                    # if subset(dominance_mask(subgraph), reachable_roots(x)) && any(reachable_variables(subgraph) .& reachable_variables(x)) #.& is 4x faster than .&& and allocates 1/12th as much. Making this one change in both versions of subgraph_exists reduces overall allocation for symbolic_jacobian! by 4x and computation time by 3x. any(a .& b) allocates, presumably to create a temporary bit vector to hold the result.

                    edge_count += 1
                end
            end

            if edge_count < 2
                return false
            end

            #     #this test appears to be wrong. Not sure why but it appears to be possible to have subgraph destroyed without having an edge at either the dominated or dominating node destroyed. 
            # dedges = filter(x -> subset(dominance_mask(subgraph), reachable_roots(x)) && any(reachable_variables(subgraph) .&& reachable_variables(x)), cedges) #need at least two constrained child edges that satisfy reachable roots and reachable variables or subgraph doesn't exist. This is a necessary but not sufficient condition for existence. Could have two properly constrained child edges and subgraph still might have been destroyed.
            # if length(dedges) < 2
            #     return false
            # elseif length(dedges) == subgraph.num_dominator_edges && length(sub_edges) == subgraph.num_dominated_edges
            #     return true #properly constrained edges for both dominator and dominated node still equal to original number of edges when subgraph was created. If the subgraph has been destroyed by factorization this cannot be true. One or the other would have to be smaller than the original value. subgraph still exists and is factorable.
            # end
        end
    end

    reclaim_edge_vector(sub_edges)

    return valid_paths(constraint, subgraph)
end

"""Returns true if the subgraph is still a factorable dominance subgraph, false otherwise"""
function subgraph_exists(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}) where {T}
    #Do fast tests that guarantee subgraph has been destroyed by factorization: no edges connected to dominated node, dominated_node or dominator node has < 2 subgraph edges
    #This is inefficient since many tests require just the number of edges but this code creates temp arrays containing the edges and then measures the length. Optimize later by having separate children and parents fields in edges structure of RnToRmGraph. Then num_parents and num_children become fast and allocation free.
    constraint = next_edge_constraint(subgraph)
    dgraph = graph(subgraph)

    g_edges = edges(dgraph)
    if get(g_edges, dominated_node(subgraph), nothing) === nothing
        return false
    else
        sub_edges = get_edge_vector()
        relation_edges!(constraint, dominated_node(subgraph), sub_edges)
        if length(sub_edges) < 2
            return false
        else
            pedges = parent_edges(dgraph, dominating_node(subgraph))
            if length(pedges) < 2
                return false
            else
                edge_count = 0

                for x in pedges
                    if subset(dominance_mask(subgraph), reachable_variables(x)) && overlap(reachable_roots(subgraph), reachable_roots(x))
                        # if subset(dominance_mask(subgraph), reachable_variables(x)) && any(reachable_roots(subgraph) .& reachable_roots(x)) #.& is 4x faster than .&& and allocates 1/12th as much. Making this one change in both versions of subgraph_exists reduces overall allocation for symbolic_jacobian! by 4x and computation time by 3x. any(a .& b) allocates, presumably to create a temporary bit vector to hold the result.
                        edge_count += 1
                    end
                end

                if edge_count < 2
                    return false
                end

                #     #this test appears to be wrong. Not sure why but it appears to be possible to have subgraph destroyed without having an edge at either the dominated or dominating node destroyed. 
                #     #     if length(dedges) == subgraph.num_dominator_edges && length(sub_edges) == subgraph.num_dominated_edges
                #     #         return true #properly constrained edges for both dominator and dominated node still equal to original number of edges when subgraph was created. subgraph still exists and is factorable.
                #     #     end
            end
        end
    end

    reclaim_edge_vector(sub_edges)

    return valid_paths(constraint, subgraph)
end


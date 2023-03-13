
# next_edge_constraint(subgraph::FactorableSubgraph{T,S}) where {T,S<:Union{PostDominatorSubgraph,DominatorSubgraph}} = next_edge_constraint(subgraph, dominance_mask(subgraph))
# next_edge_constraint(sub::FactorableSubgraph{T,DominatorSubgraph}, roots_mask::BitVector) where {T} = PathConstraint(graph(sub), true, roots_mask)
# next_edge_constraint(sub::FactorableSubgraph{T,PostDominatorSubgraph}, variables_mask::BitVector) where {T} = PathConstraint(graph(sub), false, variables_mask)

next_edge_constraint(sub::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = PathConstraint(dominating_node(sub), graph(sub), false, reachable_roots(sub), dominance_mask(sub))
next_edge_constraint(sub::FactorableSubgraph{T,DominatorSubgraph}) where {T} = PathConstraint(dominating_node(sub), graph(sub), true, dominance_mask(sub), reachable_variables(sub))
top_down_constraint(sub::FactorableSubgraph{T,DominatorSubgraph}) where {T} = PathConstraint()

"""Evaluates the subgraph, creates a new edge with this value, and then inserts the new edge into `graph`"""
function add_edge!(graph::RnToRmGraph, subgraph::FactorableSubgraph, subgraph_value::Node)
    verts = subgraph_vertices(subgraph)
    edge = PathEdge(verts[1], verts[2], subgraph_value, reachable_variables(subgraph), reachable_roots(subgraph))
    add_edge!(graph, edge)
end


function format_string(rv_string)
    tmp = ""
    for (i, rstr) in pairs(rv_string)
        if rstr != ""
            if i == lastindex(rv_string)
                tmp *= rstr
            else
                tmp *= "$rstr, "
            end
        end
    end
    if rv_string[end] == ""
        tmp = tmp[1:end-2]
    end
    return tmp
end

function print_subgraph(io, a, rv_string, subgraph_type)
    print(io, "$subgraph_type([")
    tmp = format_string(rv_string)

    print(io, tmp)
    print(io, "]")
    print(io, ", $(a.subgraph) , $(times_used(a)))")
end

# function Base.show(io::IO, a::FactorableSubgraph{T,DominatorSubgraph}) where {T}
#     root_string = [bvar == 0 ? "" : "r$i" for (i, bvar) in pairs(dom_mask(a))]
#     print_subgraph(io, a, root_string, "D")

# end

# function Base.show(io::IO, a::FactorableSubgraph{T,PostDominatorSubgraph}) where {T}
#     var_string = [bvar == 0 ? "" : "v$i" for (i, bvar) in pairs(pdom_mask(a))]
#     print_subgraph(io, a, var_string, "P")
# end

"""Holds information for factorable subgraph that is both a dom and pdom."""



"""
Finds the factorable dom subgraph associated with `node_index` if there is one.  

Given vertex `node_index` and the idom table for root `rᵢ`
find `a=idom(b)`. `a` is guaranteed to be in in the subset of the ℝⁿ->ℝᵐ graph reachable from root `rᵢ` because otherwise it would not be in the idom table associated with `rᵢ`.

The subgraph `(a,b)` is factorable if there are two or more parent edges of `b` which are both on the path to root `ri`.  

Proof: 
For `(a,b)` to be factorable there must be at least two paths from `a` downward to `b` and at least two paths upward from `b` to `a`. Because `a=idom(b)` we only need to check the latter case. If this holds then there must also be two paths downward from `a` to `b`. 

Proof:

Assume `a` has only one child. Since `a=idom(b)` all upward paths from `b` must pass through `a`. Since `a` has only one child, `nᵢ`, all paths from `b` must first pass through `nᵢ` before passing through `a`. But then `nᵢ` would be the idom of `b`, not `a`, which violates `a=idom(b)`. Hence there must be two downward paths from `a` to `b`.
"""
function dom_subgraph(graph::RnToRmGraph, root_index::Integer, dominated_node::Integer, idom)
    tmp = _node_edges(edges(graph), dominated_node)

    if tmp === nothing #no edges upward so must be a root. dominated_node can't be part of a factorable dom subgraph.
        return nothing
    else
        dominated_edges = parents(tmp)
        count = 0
        for edge in dominated_edges
            if bott_vertex(edge) == dominated_node && reachable_roots(edge)[root_index] #there is an edge upward which has a path to the root
                count += 1
                if count > 1
                    return (idom[dominated_node], dominated_node)
                end
            end
        end
        return nothing
    end
end

function pdom_subgraph(graph::RnToRmGraph, variable_index::Integer, dominated_node::Integer, pidom)
    tmp = _node_edges(edges(graph), dominated_node)
    if tmp === nothing #no edges upward so must be a root. dominated_node can't be part of a factorable pdom subgraph.
        return nothing
    else
        dominated_edges = children(tmp)

        count = 0
        for edge in dominated_edges
            if top_vertex(edge) == dominated_node && reachable_variables(edge)[variable_index] #there is an edge downward which has a path to the variable
                count += 1
                if count > 1
                    return (pidom[dominated_node], dominated_node)
                end
            end
        end
        return nothing
    end
end

#No longer used but may be useful for something in the future. TODO: delete if not used soon.
# """Necessary but not sufficient conditions for gr2 ⊂ gr1. Used to determine order to process subgraphs. Contained subgraphs must always be processed first, then sort by number of uses."""
# function maybe_contains(gr1::FactorableSubgraph, gr2::FactorableSubgraph)
#     (a, b) = gr2.subgraph
#     (c, d) = gr1.subgraph
#     if isdominator(gr2) && isdominator(gr1) #dom,dom case
#         return d < a ≤ c && d < b < c
#     elseif isdominator(gr2) && ispostdominator(gr1) #dom,pdom case
#         return c < a ≤ d && c ≤ b < d
#     elseif ispostdominator(gr2) && isdominator(gr1) #pdom,dom case
#         return d ≤ a < c && d < b ≤ c
#     else #pdom,pdom case
#         return c ≤ a < d && c < b < d
#     end
# end

"""returns true if a should be sorted before b"""
function factor_order(a::FactorableSubgraph, b::FactorableSubgraph)
    if times_used(a) > times_used(b) #num_uses of contained subgraphs always ≥ num_uses of containing subgraphs. Contained subgraphs should always be factored first. It might be that a ⊄ b, but it's still correct to factor a before b.
        return true
    elseif times_used(b) > times_used(a)
        return false
    else # if a ⊂ b then diff(a) < diff(b) where diff(x) = abs(dominating_node(a) - dominated_node(a)). Might be that a ⊄ b but it's safe to factor a first and if a ⊂ b then it will always be factored first.
        diffa = node_difference(a)
        diffb = node_difference(b)

        if diffa < diffb
            return true
        else
            return false #can factor a,b in either order 
        end
    end
end
export factor_order

sort_in_factor_order!(a::AbstractVector{T}) where {T<:FactorableSubgraph} = sort!(a, lt=factor_order)
export sort_in_factor_order!

"""
Given subgraph `(a,b)` in the subset of the ℝⁿ->ℝᵐ graph reachable from root `rᵢ`.  
`(a,b)` is factorable iff: 

`a > b`  
`&& dom(a,b) == true`  
`&& num_parents(b) > 1 for parents on the path to root `rᵢ` through `a`  
`&& num_children(a) > 1 for children on the path to `b`  
 
or

subgraph `(a,b)` is in the subset of the ℝⁿ->ℝᵐ graph reachable from variable `vⱼ` 

`a < b`  
`&& pdom(a,b) == true`  
`&& num_parents(a) > 1 for parents on the path to `b`  
`&& num_children(b) > 1 for children on the path to variable `vᵢ` through `a`  
"""
function compute_factorable_subgraphs(graph::RnToRmGraph{T}) where {T}
    dom_subgraphs = Dict{Tuple{T,T},BitVector}()
    pdom_subgraphs = Dict{Tuple{T,T},BitVector}()
    dual_subgraphs = Tuple{T,T,BitVector,BitVector}[]
    idoms = compute_dominance_tables(graph, true)
    pidoms = compute_dominance_tables(graph, false)
    subgraph_map = Dict{Tuple{T,T},Tuple{FactorableSubgraph,T}}()

    for root_index in 1:codomain_dimension(graph)
        for dominated in keys(idoms[root_index])
            dsubgraph = dom_subgraph(graph, root_index, dominated, idoms[root_index])
            if dsubgraph !== nothing
                existing_subgraph = get(dom_subgraphs, dsubgraph, nothing)
                if existing_subgraph !== nothing
                    dom_subgraphs[dsubgraph][root_index] = 1
                else
                    tmp = falses(codomain_dimension(graph))
                    tmp[root_index] = 1
                    dom_subgraphs[dsubgraph] = tmp
                end
            end
        end
    end

    for variable_index in 1:domain_dimension(graph)
        for dominated in keys(pidoms[variable_index])
            psubgraph = pdom_subgraph(graph, variable_index, dominated, pidoms[variable_index])
            if psubgraph !== nothing
                existing_p_subgraph = get(pdom_subgraphs, psubgraph, nothing)

                if existing_p_subgraph !== nothing
                    pdom_subgraphs[psubgraph][variable_index] = 1
                else
                    tmp = falses(domain_dimension(graph))
                    tmp[variable_index] = 1
                    pdom_subgraphs[psubgraph] = tmp
                end
            end
        end
    end


    #convert to factorable subgraphs

    result = Vector{FactorableSubgraph}(undef, length(dom_subgraphs) + length(pdom_subgraphs)) #slightly faster and more memory efficient than pushing to a zero length array
    empty!(result)

    #Explanation of the computation of uses. Assume key[1] > key[2] so subgraph is a dom. subgraphs[key] stores the number of roots for which this dom was found to be factorable. 
    #For each root that has the dom as a factorable subgraph the number of paths from the bottom node of the subgraph to the variables will be the same. Total number of uses is the
    #product of number of roots with dom factorable * number of paths to variables. Similar argument holds for pdom factorable subgraphs.
    for key in keys(dom_subgraphs)
        dominator = key[1]
        dominated = key[2]
        subgraph = dominator_subgraph(graph, dominator, dominated, dom_subgraphs[key], reachable_roots(graph, dominator), reachable_variables(graph, dominated))

        push!(result, subgraph)
        subgraph_map[(dominator, dominated)] = (subgraph, lastindex(result))
    end

    for key in keys(pdom_subgraphs)
        dominator = key[1]
        dominated = key[2]
        subgraph = postdominator_subgraph(graph, dominator, dominated, pdom_subgraphs[key], reachable_roots(graph, dominated), reachable_variables(graph, dominator))

        push!(result, subgraph)
        subgraph_map[(dominator, dominated)] = (subgraph, lastindex(result))
    end

    sort_in_factor_order!(result) #return subgraphs sorted in lexicographic order. If the subgraphs are processed from lowest index to highest index this will guarantee that innermost graphs are processed first.

    return result, subgraph_map
end
export compute_factorable_subgraphs

function multiply_sequence(path::AbstractVector{S}) where {S<:PathEdge}
    if length(path) == 1
        return edge_value(path[1])
    end

    run = Node[]
    count = 2
    prod = Node(1.0)
    run_start = times_used(path[1])
    push!(run, edge_value(path[1]))

    for val in @view path[2:end]
        if times_used(val) != run_start
            runprod = Node(1.0)
            for runval in run
                runprod *= runval
            end
            empty!(run)
            prod *= runprod
            run_start = times_used(val)
        end

        push!(run, edge_value(val))

        if count == length(path)
            runprod = Node(1.0)
            for runval in run
                runprod *= runval
            end
            prod *= runprod
        end
        count += 1
    end
    return prod
end
export multiply_sequence


function path_sort_order(x, y)
    if times_used(x) > times_used(y)
        return true
    elseif times_used(x) < times_used(y)
        return false
    else
        return top_vertex(x) > top_vertex(y)
    end
end
export path_sort_order

const EDGE_CACHE = Vector{PathEdge{Int64}}[]
peak_cache_size = 0

function get_edge_vector()
    if length(EDGE_CACHE) != 0
        tmp = pop!(EDGE_CACHE)
        empty!(tmp)
    else
        PathEdge{Int64}[]
    end
end

function reclaim_edge_vector(edges::Vector{PathEdge{Int64}})
    global peak_cache_size
    push!(EDGE_CACHE, edges)
    if length(EDGE_CACHE) > peak_cache_size
        peak_cache_size = length(EDGE_CACHE)
    end
    return nothing
end

function evaluate_subgraph(subgraph::FactorableSubgraph{T,S}) where {T,S<:Union{DominatorSubgraph,PostDominatorSubgraph}}
    constraint = next_edge_constraint(subgraph)
    sum = Node(0.0)

    rel_edges = get_edge_vector()
    for edge in relation_edges(constraint, dominated_node(subgraph), rel_edges)
        pedges = get_edge_vector()
        pedges = edges_on_path(constraint, dominating_node(subgraph), S == DominatorSubgraph, edge)
        #sort by num_uses then from largest to smallest postorder number
        @assert pedges !== nothing

        sort!(pedges, lt=path_sort_order)
        sum += multiply_sequence(pedges)
        #find sequences of equal times_used and multiply them. Then multiply each of the collapsed sequences to get the final result
        reclaim_edge_vector(pedges)
    end
    reclaim_edge_vector(rel_edges)
    return sum
end
export evaluate_subgraph

make_factored_edge(subgraph::FactorableSubgraph{T,DominatorSubgraph}) where {T} = PathEdge(dominating_node(subgraph), dominated_node(subgraph), evaluate_subgraph(subgraph), reachable_variables(subgraph), dominance_mask(subgraph))
export make_factored_edge
make_factored_edge(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = PathEdge(dominating_node(subgraph), dominated_node(subgraph), evaluate_subgraph(subgraph), dominance_mask(subgraph), reachable_roots(subgraph))
export make_factored_edge

struct MaskableEdge{S<:Union{DominatorSubgraph,PostDominatorSubgraph}}
    edge::PathEdge
    mask::BitVector
end

function reset_reachable(me::MaskableEdge{DominatorSubgraph})
    zero_roots!(me.edge, me.mask)
    if is_zero(reachable_roots(me.edge))
        return me.edge
    else
        return nothing
    end
end

function reset_reachable(me::MaskableEdge{PostDominatorSubgraph})
    zero_variables!(me.edge, me.mask)
    if is_zero(reachable_variables(me.edge))
        return me.edge
    else
        return nothing
    end
end


#not efficient. Make code work correctly then optimize.
function edges_on_path(next_node_constraint, dominating::T, is_dominator::Bool, current_edge, result::Union{Nothing,Vector{PathEdge{Int64}}}=nothing) where {T}
    if result === nothing
        result = get_edge_vector()
        empty!(result)
    else
        result = PathEdge{T}[]
    end

    while true
        push!(result, current_edge)
        if is_dominator && top_vertex(current_edge) == dominating
            break
        end
        if !is_dominator && bott_vertex(current_edge) == dominating
            break
        end
        tmp = get_edge_vector()
        edge_array = relation_edges(next_node_constraint, current_edge, tmp)

        if length(edge_array) == 0 || length(edge_array) ≥ 2
            return nothing #there is either a branch in the edge path or there is no edge beyond current_edge that leads to the dominating node. This can only occur if the subgraph has been destroyed by factorization so stop processing
        end


        current_edge = edge_array[1]
        reclaim_edge_vector(edge_array)
    end
    return result
end
export edges_on_path

function compute_Vset(constraint::PathConstraint{T}, dominating_node::T, dominated_node::T) where {T}
    Vset = trues(domain_dimension(graph(constraint)))
    tmp = get_edge_vector()
    relation_edges(constraint, dominated_node, tmp)
    for start_edge in tmp
        pedges = edges_on_path(constraint, dominating_node, true, start_edge)

        #reset roots in V, if possible. All edges higher in the path than the first vertex with more than one child cannot be reset.
        for pedge in pedges
            Vset .= Vset .& reachable_variables(pedge)
        end
    end
    reclaim_edge_vector(tmp)
    return Vset
end

"""caller has to pass in graph and edges_to_delete vector"""
function factor_subgraph!(subgraph::FactorableSubgraph{T,S}, sub_eval::Union{Nothing,Node}) where {T,S<:DominatorSubgraph}
    a = graph(subgraph)
    @assert subgraph_exists(subgraph)
    Rdom = dominance_mask(subgraph) #roots for which subgraph is factorable, i.e., dominator(subgraph) dominates every vertex in the subgraph.
    edges_to_reset = MaskableEdge{S}[]
    Vset = trues(domain_dimension(graph(subgraph)))

    R_edges = next_edge_constraint(subgraph)
    dominator = dominating_node(subgraph)
    dominated = dominated_node(subgraph)
    Vset = compute_Vset(R_edges, dominator, dominated)

    for start_edge in relation_edges(R_edges, dominated)
        pedges = edges_on_path(R_edges, dominator, true, start_edge)

        #reset roots in R, if possible. All edges higher in the path than the first vertex with more than one child cannot be reset.
        for pedge in pedges
            Vval = set_diff(reachable_variables(pedge), Vset)
            if !is_zero(Vval) #need to create and edge to accomodate the part of the reachable set that is not in Rset
                add_edge!(a, PathEdge(top_vertex(pedge), bott_vertex(pedge), edge_value(pedge), Vval, Rdoml)) #this edge only has reachable roots outside Rset. Need to add this here rather than in factor_one_subgraph because dual processing may need to look at these edges
                mask_variables!(pedge, Vset) #this edge only has the reachable roots in Rset
            end

            if roots_resettable(pedge, Rdom)
                push!(edges_to_reset, MaskableEdge{DominatorSubgraph}(pedge, Rdom))
            end
            if top_vertex(pedge) == dominator || length(child_edges(a, top_vertex(pedge))) > 1
                break
            end
        end
    end

    if sub_eval === nothing
        return edges_to_reset, nothing
    else
        return edges_to_reset, PathEdge(dominator, dominated, sub_eval, Vset, Rdom)  #return edges marked for resetting and new added edge that must be added to the graphend
    end
end
export factor_subgraph!

function compute_Rset(constraint::PathConstraint{T}, dominating_node::T, dominated_node::T) where {T}
    Rset = trues(codomain_dimension(graph(constraint)))
    tmp = get_edge_vector()
    relation_edges(constraint, dominated_node, tmp)
    for start_edge in tmp
        pedges = edges_on_path(constraint, dominating_node, false, start_edge)

        #reset roots in V, if possible. All edges higher in the path than the first vertex with more than one child cannot be reset.
        for pedge in pedges
            Rset .= Rset .& reachable_roots(pedge)
        end
    end
    reclaim_edge_vector(tmp)
    return Rset
end

"""caller has to pass in graph and edges_to_delete vector. There is much redundancy with the DominatorSubgraph version of this code. But the differences are numerous (13 different locations) and the code is easier to read if split into two functions."""
function factor_subgraph!(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}, sub_eval::Union{Nothing,Node}) where {T}
    a = graph(subgraph)
    @assert subgraph_exists(subgraph)
    Vdom = dominance_mask(subgraph) #roots for which subgraph is factorable, i.e., dominator(subgraph) dominates every vertex in the subgraph.
    edges_to_reset = MaskableEdge{PostDominatorSubgraph}[]

    V_edges_up = next_edge_constraint(subgraph)
    dominator = dominating_node(subgraph)
    dominated = dominated_node(subgraph)

    Rset = compute_Rset(V_edges_up, dominator, dominated) #these are the roots for which this subgraph is factorable, i.e., only paths that include these roots will get the factored value of the subgraph.

    for start_edge in relation_edges(V_edges_up, dominated)
        pedges = edges_on_path(V_edges_up, dominator, false, start_edge)

        #reset roots in V, if possible. All edges higher in the path than the first vertex with more than one child cannot be reset.
        for pedge in pedges
            Rval = set_diff(reachable_roots(pedge), Rset)
            if !is_zero(Rval) #need to create an edge to accomodate the part of the reachable set that is not in Rset
                add_edge!(a, PathEdge(top_vertex(pedge), bott_vertex(pedge), edge_value(pedge), Vdom, Rval)) #this edge only has reachable roots outside Rset. Need to add this here rather than in factor_one_subgraph because dual processing may need to look at these edges
                mask_roots!(pedge, Rset) #this edge only has the reachable roots in Rset
            end

            if variables_resettable(pedge, Vdom)
                push!(edges_to_reset, MaskableEdge{PostDominatorSubgraph}(pedge, Vdom))
            end
            if bott_vertex(pedge) == dominator || length(parent_edges(a, bott_vertex(pedge))) > 1 || is_root(a, bott_vertex(pedge)) #is_root special case for post dominator subgraphs since a root can be a child of another root. Variables cannot have this relationship. If the bottom vertex of the edge is a root that means there is a path to a root that does not go through dominated. All edges below this point in the graph cannot have their variable bits reset.
                break
            end
        end
    end
    if sub_eval === nothing
        return edges_to_reset, nothing
    else
        return edges_to_reset, PathEdge(dominator, dominated, sub_eval, Vdom, Rset)  #return edges marked for resetting and new added edge that must be added to the graph
    end
end


function print_edges(a, msg)
    println(msg)
    for edge in edges(a)
        println(edge)
    end
end

"""call this function when the first subgraph processed in a dual pair is a DominatorSubgraph"""
function dual_subgraph(first_graph::FactorableSubgraph{T,DominatorSubgraph}, second_graph::FactorableSubgraph{T,PostDominatorSubgraph}) where {T}
    @assert first_graph.graph === second_graph.graph

    graph = first_graph.graph
    Vdom = dominance_mask(second_graph)
    R = reachable_roots(first_graph)
    Rdom = dominance_mask(first_graph)
    R̄ = set_diff(R, Rdom)
    overlap_graph = postdominator_subgraph(graph, dominating_node(second_graph), dominated_node(second_graph), Vdom, Rdom, Vdom)
    if is_zero(R̄)
        return overlap_graph, nothing
    else
        return overlap_graph, postdominator_subgraph(graph, dominating_node(second_graph), dominated_node(second_graph), Vdom, R̄, Vdom) #could let the second Vdom be reachable_variables(second_graph) because this term is not used in computing reachable set for new factored subgraph edge. But safer to do this in case factor_subgraph! code is changed. Then guaranteed not to accidentally get too wide reachability.
    end
end

"""call this function when the first subgraph processed in a dual pair is a PostDominatorSubgraph"""
function dual_subgraph(first_graph::FactorableSubgraph{T,PostDominatorSubgraph}, second_graph::FactorableSubgraph{T,DominatorSubgraph}) where {T}
    @assert first_graph.graph === second_graph.graph

    graph = first_graph.graph
    Rdom = dominance_mask(second_graph)
    R = reachable_roots(first_graph)
    V = reachable_variables(first_graph)
    Vdom = dominance_mask(first_graph)
    V̄ = set_diff(V, Vdom)
    R̄ = set_diff(R, Rdom)
    overlap_graph = dominator_subgraph(graph, dominating_node(second_graph), dominated_node(second_graph), Rdom, Rdom, Vdom)
    if is_zero(V̄)
        return overlap_graph, nothing #complete overlap between dual graphs so no extra edges to be processed in the V̄ graph.
    else
        return overlap_graph, dominator_subgraph(graph, dominating_node(second_graph), dominated_node(second_graph), Rdom, Rdom, V̄) #have a dual residual graph that must be processed later
    end
end

function process_dual_subgraph(subgraph::FactorableSubgraph{T}, subgraph_dict, subgraph_list) where {T}
    dual_vertices = reverse(subgraph_vertices(subgraph))
    dual_info = get(subgraph_dict, dual_vertices, nothing)

    dual_edges_to_reset = nothing
    if dual_info !== nothing # there is a dual subgraph so have more processing to do
        dual_graph, list_index = dual_info
        delete!(subgraph_dict, dual_vertices) #delete original dual subgraph so don't attempt to process it when it is encountered later in the subgraph_list

        if subgraph_exists(dual_graph) #subgraph may have been destroyed so make sure it still exists before processing
            overlap_dual, residual_dual = dual_subgraph(subgraph, dual_graph) #find the part of the dual graph that must be processed now with subgraph and the residual part that will be processed later.

            if residual_dual !== nothing
                subgraph_dict[dual_vertices] = (residual_dual, list_index) #substitute new residual dual graph in place of the old one
                subgraph_list[list_index] = residual_dual
            end

            dual_edges_to_reset, _ = factor_subgraph!(overlap_dual, nothing)
            return dual_edges_to_reset
        else
            return PathEdge{T}[]
        end
    end
end


function factor_one_subgraph!(a::RnToRmGraph, subgraph::FactorableSubgraph, subgraph_list, subgraph_dict)
    if get(subgraph_dict, subgraph_vertices(subgraph), nothing) !== nothing #subgraph may already have been processed in which case it will have been removed from subgraph_dict. This happens with dual subgraphs (a,b), (b,a). When (a,b) is processed it will also process (b,a) and (b,a) may be removed from subgraph_dict if it completely overlaps with (a,b).
        # @info("before processing subgraph $subgraph")

        if subgraph_exists(subgraph)
            sub_eval = evaluate_subgraph(subgraph)
            edges_to_delete = PathEdge[]
            edges_to_reset, factored_edge_to_add = factor_subgraph!(subgraph, sub_eval)
            dual_edges_to_reset = process_dual_subgraph(subgraph, subgraph_dict, subgraph_list)
            combined = dual_edges_to_reset !== nothing ? vcat(edges_to_reset, dual_edges_to_reset) : edges_to_reset

            for edge in combined
                tmp = reset_reachable(edge)
                if tmp !== nothing
                    push!(edges_to_delete, tmp)
                end
            end

            add_edge!(a, factored_edge_to_add)

            delete_edge!.(Ref(a), unique(edges_to_delete)) #dom,pdom may both add an edge to edges_to_delete. This will cause an assertion failure so make sure only delete edges once.
        end

        delete!(subgraph_dict, subgraph) #not strictly necessary because should never encounter the same subgraph twice so could leave the subgraph in the dictionary. But makes debugging easier when looking at small test cases.
    end

    return nothing #return nothing so people don't mistakenly think this is returning a copy of the original graph
end

function factor!(a::RnToRmGraph{T}) where {T}
    subgraph_list, subgraph_dict = compute_factorable_subgraphs(a)

    for subgraph in subgraph_list
        factor_one_subgraph!(a, subgraph, subgraph_list, subgraph_dict)
    end

    return nothing #return nothing so people don't mistakenly think this is returning a copy of the original graph
end
export factor!

function follow_path(a::RnToRmGraph{T}, root_index::Integer, var_index::Integer) where {T}
    current_node_index = root_index_to_postorder_number(a, root_index)
    path_product = PathEdge{T}[]

    while true
        curr_edges = filter(x -> is_root_reachable(x, root_index) && is_variable_reachable(x, var_index), child_edges(a, current_node_index))
        if length(curr_edges) == 0
            break
        else
            @assert length(curr_edges) == 1 "Should only be one path from root $root_index to variable $var_index. Instead have $(length(curr_edges)) children from a node on the path"
            push!(path_product, curr_edges[1])
            current_node_index = bott_vertex(curr_edges[1])
        end
    end
    if length(path_product) == 0
        product = Node(0.0)
    else
        sort!(path_product, lt=((x, y) -> num_uses(x) > num_uses(y))) #sort larger num uses edges first
        product = Node(1.0)
        for term in path_product
            product *= edge_value(term)
        end
    end
    return product
end

function evaluate_path(graph::RnToRmGraph, root_index::Integer, var_index::Integer)
    node_value = root(graph, root_index)
    if !is_tree(node_value) #root contains a variable or constant
        if is_variable(node_value)
            if variable(graph, var_index) == node_value
                return 1.0 #taking a derivative with respect to itself, which is 1. Need to figure out a better way to get the number type right
            else
                return 0.0 #taking a derivative with respect to a different variable, which is 0.
            end
        else
            return 0.0 #root is a constant
        end
    else #root contains a graph which has been factored so that there should be a single linear path from each root to each variable with no branching
        return follow_path(graph, root_index, var_index)
    end
end

"""deletes edges which do not have a path to a variable"""
function path_to_variable!(graph, current_edge)
    if is_variable(graph, bott_vertex(current_edge))
        return true
    elseif is_constant(graph, bott_vertex(current_edge))
        delete_edge!(graph, current_edge, true)
        return false
    else
        is_path = false

        edge_copy = copy(child_edges(graph, current_edge))
        for cedge in edge_copy #have to use a copy because the graph data structure referenced by child_edges will be mutated
            tmp = path_to_variable!(graph, cedge)
            is_path = is_path || tmp
        end

        if is_path == false
            delete_edge!(graph, current_edge, true)
            return false
        end
        return is_path
    end
end

function remove_dangling_edges(graph::RnToRmGraph)
    #might be legal for root to have multiple dangling paths but I don't think any other nodes should. Requires proof, might not be true.
    for root_index in 1:codomain_dimension(graph)
        edge_copy = copy(child_edges(graph, root_index_to_postorder_number(graph, root_index)))
        for cedge in edge_copy #need to copy edge data structure because the child edge array in the graph is being edited by path_to_variable!. This causes premature termination of the loop.
            path_to_variable!(graph, cedge)
        end
    end
end


"""Factors the graph then computes jacobian matrix. Destructive."""
function symbolic_jacobian!(a::RnToRmGraph, variable_ordering::AbstractVector{T}) where {T<:Node}
    indim = domain_dimension(a)
    outdim = codomain_dimension(a)

    result = Matrix{Node}(undef, outdim, indim)
    factor!(a)

    remove_dangling_edges(a)

    for (i, var) in pairs(variable_ordering)
        var_index = variable_node_to_index(a, var)
        for root_index in 1:codomain_dimension(a)
            result[root_index, i] = evaluate_path(a, root_index, var_index)
        end
    end

    return result
end
export symbolic_jacobian!

symbolic_jacobian!(a::RnToRmGraph) = symbolic_jacobian!(a, variables(a))


jacobian_function!(graph::RnToRmGraph) = jacobian_function!(graph, variables(graph))

function jacobian_function!(graph::RnToRmGraph, variable_order::AbstractVector{S}) where {S<:Node}
    tmp = symbolic_jacobian!(graph, variable_order)

    node_to_var = Dict{Node,Union{Symbol,Real}}()
    body = Expr(:block)
    push!(body.args, :(result = fill(0.0, $(size(tmp)))))
    all_vars = variables(graph)
    if variable_order === nothing
        ordering = all_vars
    else
        ordering = Node.(variable_order)
    end
    @assert Set(all_vars) ⊆ Set(ordering) "Not every variable in the graph had a corresponding ordering variable."

    for (i, node) in pairs(tmp)
        node_body, variable = function_body(node, node_to_var)
        push!(node_body.args, :(result[$i] = $variable))
        push!(body.args, node_body)
    end

    push!(body.args, :(return result))

    return @RuntimeGeneratedFunction(Expr(:->, Expr(:tuple, map(x -> node_symbol(x), ordering)...), body))
end
export jacobian_function!

"""Returns the number of unique nodes in a jacobian. Used to roughly estimate number of operations to evaluation jacobian"""
function num_unique_nodes(jacobian::Matrix{Node})
    nodes = Set{Node}()
    for oned in all_nodes.(jacobian)
        union!(nodes, oned)
    end
    return length(nodes)
end
export num_unique_nodes

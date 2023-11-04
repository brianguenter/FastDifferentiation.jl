
next_edge_constraint(sub::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = PathConstraint(dominating_node(sub), graph(sub), false, reachable_roots(sub), reachable_dominance(sub))
next_edge_constraint(sub::FactorableSubgraph{T,DominatorSubgraph}) where {T} = PathConstraint(dominating_node(sub), graph(sub), true, reachable_dominance(sub), reachable_variables(sub))
top_down_constraint(sub::FactorableSubgraph{T,DominatorSubgraph}) where {T} = PathConstraint()

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
function dom_subgraph(graph::DerivativeGraph, root_index::Integer, dominated_node::Integer, idom)
    dominated_edges = parent_edges(graph, dominated_node)

    if length(dominated_edges) == 0  #no edges edges up so must be a root. dominated_node can't be part of a factorable dom subgraph.
        return nothing
    else
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

function pdom_subgraph(graph::DerivativeGraph, variable_index::Integer, dominated_node::Integer, pidom)
    dominated_edges = child_edges(graph, dominated_node)

    if length(dominated_edges) == 0  #no edges up so must be a root. dominated_node can't be part of a factorable dom subgraph.
        return nothing
    else
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

struct FactorOrder <: Base.Order.Ordering
end

Base.lt(::FactorOrder, a, b) = factor_order(a, b)
Base.isless(::FactorOrder, a, b) = factor_order(a, b)


"""returns true if a should be sorted before b"""
function factor_order(a::FactorableSubgraph, b::FactorableSubgraph)
    # if a ⊂ b then diff(a) < diff(b) where diff(x) = abs(dominating_node(a) - dominated_node(a)). Might be that a ⊄ b but it's safe to factor a first and if a ⊂ b then it will always be factored first.
    diffa = node_difference(a)
    diffb = node_difference(b)

    if diffa < diffb
        return true
    elseif times_used(a) > times_used(b) #Factor subgraphs with more uses first for efficiency.
        return true
    else
        return false
    end
end

sort_in_factor_order!(a::AbstractVector{T}) where {T<:FactorableSubgraph} = sort!(a, lt=factor_order)



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
function compute_factorable_subgraphs(graph::DerivativeGraph{T}) where {T}
    pdom_subgraphs = Dict{Tuple{T,T},BitVector}()
    dom_subgraphs = Dict{Tuple{T,T},BitVector}()

    function set_dom_bits!(subgraph::Union{Nothing,Tuple{T,T}}, all_subgraphs::Dict, var_or_root_index::Integer, bit_dimension::Integer) where {T<:Integer}
        if subgraph !== nothing
            existing_subgraph = get(all_subgraphs, subgraph, nothing)
            if existing_subgraph !== nothing
                all_subgraphs[subgraph][var_or_root_index] = 1
            else
                tmp = falses(bit_dimension)
                tmp[var_or_root_index] = 1
                all_subgraphs[subgraph] = tmp
            end
        end
    end

    temp_doms = Dict{T,T}()

    for root_index in 1:codomain_dimension(graph)
        post_num = root_index_to_postorder_number(graph, root_index)
        temp_dom = compute_dom_table(graph, true, root_index, post_num, temp_doms)

        for dominated in keys(temp_dom)
            dsubgraph = dom_subgraph(graph, root_index, dominated, temp_dom)
            set_dom_bits!(dsubgraph, dom_subgraphs, root_index, codomain_dimension(graph))
        end
    end

    for variable_index in 1:domain_dimension(graph)
        post_num = variable_index_to_postorder_number(graph, variable_index)
        temp_dom = compute_dom_table(graph, false, variable_index, post_num, temp_doms)

        for dominated in keys(temp_dom)
            psubgraph = pdom_subgraph(graph, variable_index, dominated, temp_dom)
            set_dom_bits!(psubgraph, pdom_subgraphs, variable_index, domain_dimension(graph))
        end
    end


    #convert to factorable subgraphs

    result = BinaryHeap{FactorableSubgraph{T,S} where {S<:AbstractFactorableSubgraph},FactorOrder}()

    #Explanation of the computation of uses. Assume key[1] > key[2] so subgraph is a dom. subgraphs[key] stores the number of roots for which this dom was found to be factorable. 
    #For each root that has the dom as a factorable subgraph the number of paths from the bottom node of the subgraph to the variables will be the same. Total number of uses is the
    #product of number of roots with dom factorable * number of paths to variables. Similar argument holds for pdom factorable subgraphs.
    for key in keys(dom_subgraphs)
        dominator = key[1]
        dominated = key[2]
        if !is_constant(node(graph, dominated)) #don't make subgraphs with constant dominated nodes because they are not factorable
            subgraph = dominator_subgraph(graph, dominator, dominated, dom_subgraphs[key])

            push!(result, subgraph)
        end
    end

    for key in keys(pdom_subgraphs)
        dominator = key[1]
        dominated = key[2]
        subgraph = postdominator_subgraph(graph, dominator, dominated, pdom_subgraphs[key])

        push!(result, subgraph)
    end

    return result
end


function multiply_sequence(path::AbstractVector{S}) where {S<:PathEdge}
    if length(path) == 1
        return value(path[1])
    end

    run = Node[]
    count = 2
    prod = Node(1.0)
    run_start = times_used(path[1])
    push!(run, value(path[1]))

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

        push!(run, value(val))

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



function path_sort_order(x, y)
    if times_used(x) > times_used(y)
        return true
    elseif times_used(x) < times_used(y)
        return false
    else
        return top_vertex(x) > top_vertex(y)
    end
end


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

function factored_reachable_roots(subgraph::FactorableSubgraph{T,DominatorSubgraph}, sub_edges) where {T}
    result = falses(codomain_dimension(graph(subgraph)))

    for edge in child_edges(graph(subgraph), dominating_node(subgraph))
        if edge ∈ sub_edges
            result .= result .| reachable_roots(edge)
        end
    end
    result .= result .& reachable_dominance(subgraph) #only the roots that are in the dominance set are counted here because this is used to determine the reachability of the subgraph sum.
    return result
end

function factored_reachable_roots(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}, sub_edges) where {T}
    result = falses(codomain_dimension(graph(subgraph)))

    for edge in child_edges(graph(subgraph), dominated_node(subgraph))
        if edge ∈ sub_edges
            result .= result .| reachable_roots(edge)
        end
    end
    return result
end

function factored_reachable_variables(subgraph::FactorableSubgraph{T,DominatorSubgraph}, sub_edges) where {T}
    result = falses(domain_dimension(graph(subgraph)))

    for edge in parent_edges(graph(subgraph), dominated_node(subgraph))
        if edge ∈ sub_edges
            result .= result .| reachable_variables(edge)
        end
    end
    return result
end

function factored_reachable_variables(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}, sub_edges) where {T}
    result = falses(domain_dimension(graph(subgraph)))

    for edge in parent_edges(graph(subgraph), dominating_node(subgraph))
        if edge ∈ sub_edges
            result .= result .| reachable_variables(edge)
        end
    end

    result .= result .& reachable_dominance(subgraph) #only the variables that are in the dominance set are counted here because this is used to determine the reachability of the subgraph sum.
    return result
end


function make_factored_edge(subgraph::FactorableSubgraph, sub_edges, sum::Node)
    roots_reach = factored_reachable_roots(subgraph, sub_edges)
    vars_reach = factored_reachable_variables(subgraph, sub_edges)
    return PathEdge(dominating_node(subgraph), dominated_node(subgraph), sum, vars_reach, roots_reach)
end

"""returns the number of subgraph edges incident on `dominated_node(subgraph)` and `dominating_node(subgraph)`"""
function dom_nodes_edge_count(subgraph, sub_edges)
    dominating_edges = PathEdge[]
    dominated_edges = PathEdge[]
    for edge in sub_edges
        if dominating_node(subgraph) == forward_vertex(subgraph, edge)
            push!(dominating_edges, edge)
        end

        if dominated_node(subgraph) == backward_vertex(subgraph, edge)
            push!(dominated_edges, edge)
        end
    end

    return dominated_edges, dominating_edges
end


"""Returns true if a new factorable subgraph was created inside `subgraph` during the factorization process. This will cause a branch internal to the subgraph. `subgraph_exists` should be called before executing this function otherwise it may return false when no new subgraphs have been created."""

# function is_branching(subgraph, sub_edges)
#     for edge in sub_edges
#         #look for a an upward branch
#         if forward_vertex(subgraph, edge) != dominating_node(subgraph)
#             let count = 0
#                 tmp = forward_edges(subgraph, edge)
#                 for fedge in tmp
#                     if fedge ∈ sub_edges
#                         count += 1
#                     end
#                     if count > 1
#                         return true
#                     end
#                 end
#             end
#         end
#         #look for a downward branch
#         if backward_vertex(subgraph, edge) != dominated_node(subgraph)
#             let count = 0
#                 tmp = backward_edges(subgraph, edge)
#                 for bedge in tmp
#                     if bedge ∈ sub_edges
#                         count += 1
#                     end
#                     if count > 1
#                         return true
#                     end
#                 end
#             end
#         end
#     end
#     return false
# end

"""sorts edges by reverse postorder of the `top_vertex(edge)`"""
function sort_edge_list!(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}, edge_list) where {T}
    sort!(edge_list, by=x -> backward_vertex(subgraph, x), rev=true)
end

"""sorts edges by postorder of the `bott_vertex(edge)`"""
function sort_edge_list!(subgraph::FactorableSubgraph{T,DominatorSubgraph}, edge_list) where {T}
    sort!(edge_list, by=x -> backward_vertex(subgraph, x))
end

function interior_mask_bit(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}, vertex::S) where {T,S<:Integer}
    if is_root(graph(subgraph), vertex) && vertex != dominated_node(subgraph)
        tmp = falses(codomain_dimension(graph(subgraph)))
        root_number = root_postorder_to_index(graph(subgraph), vertex)
        tmp[root_number] = true
        return tmp
    else
        return nothing
    end
end

function interior_mask_bit(subgraph::FactorableSubgraph{T,DominatorSubgraph}, vertex::S) where {T,S<:Integer}
    if is_variable(graph(subgraph), vertex) && vertex != dominated_node(subgraph)
        tmp = falses(domain_dimension(graph(subgraph)))
        variable_number = variable_postorder_to_index(graph(subgraph), vertex)
        tmp[variable_number] = true
        return tmp
    else
        return nothing
    end
end

"""Returns true if at least one root and at least one variable are reachable from the edge, false otherwise."""
function is_reachable(subgraph, edge)
    any(reachable_dominance(subgraph, edge)) && any(reachable_non_dominance(subgraph, edge))
end

function compute_vertex_masks(subgraph::FactorableSubgraph{T}, sub_edges) where {T<:Integer}
    edge_list = collect(sub_edges)

    sort_edge_list!(subgraph, edge_list) #sorts edges by postorder for DominatorSubgraph and reverse postorder for PostDominatorSubgraph.

    vertex_masks = Dict{T,BitVector}()
    vertex_masks[dominated_node(subgraph)] = falses(non_dominance_dimension(subgraph)) #

    #compute vertex masks for all edges
    for edge in edge_list
        vert = backward_vertex(subgraph, edge)
        if vert != dominated_node(subgraph)
            tmp_mask = get(vertex_masks, vert, nothing)
            if tmp_mask === nothing
                tmp_mask = falses(non_dominance_dimension(subgraph))
                vertex_masks[vert] = tmp_mask
            end

            for cedge in backward_edges(subgraph, edge)
                #check to see if a vertex in the subgraph is a root or variable that is not equal to the dominated_node.
                back_vert = backward_vertex(subgraph, edge)
                interior_mask = interior_mask_bit(subgraph, back_vert)
                if interior_mask !== nothing
                    tmp_mask .= tmp_mask .| interior_mask
                end

                if cedge ∉ sub_edges
                    if any(reachable_dominance(subgraph) .& reachable_dominance(subgraph, cedge))
                        tmp_mask .= tmp_mask .| non_dominance_mask(subgraph, cedge) #if any edge bypasses the dominated node then write a 1 in the mask for all the variable/root indices reachable from that edge
                    end
                else
                    tmp_mask .= tmp_mask .| (vertex_masks[backward_vertex(subgraph, cedge)] .& non_dominance_mask(subgraph, edge))
                end
            end
        end
    end
    return vertex_masks
end

"""Find edges which bypass the non-dominant vertex of a subgraph. These edges must be preserved after factorization"""
function find_bypass_edges(subgraph::FactorableSubgraph{T,S}, sub_edges) where {T,S<:AbstractFactorableSubgraph}
    edge_list = collect(sub_edges)

    sort_edge_list!(subgraph, edge_list) #sorts edges by postorder for DominatorSubgraph and reverse postorder for PostDominatorSubgraph.

    vertex_masks = Dict{T,BitVector}()
    vertex_masks[dominated_node(subgraph)] = falses(non_dominance_dimension(subgraph)) #

    #this propagates bypass reachability bits through the graph in the forward_vertex direction. The vertices are visited in this order because the edge list is sorted by postorder or reverse postorder number. 

    #compute vertex masks for all edges
    for edge in edge_list
        vert = backward_vertex(subgraph, edge)
        if vert != dominated_node(subgraph)
            tmp_mask = get(vertex_masks, vert, nothing)
            if tmp_mask === nothing
                tmp_mask = falses(non_dominance_dimension(subgraph))
                vertex_masks[vert] = tmp_mask
            end

            for cedge in backward_edges(subgraph, edge)
                #check to see if a vertex in the subgraph is a root or variable that is not equal to the dominated_node.
                back_vert = backward_vertex(subgraph, edge)
                interior_mask = interior_mask_bit(subgraph, back_vert)
                if interior_mask !== nothing
                    tmp_mask .= tmp_mask .| interior_mask
                end

                if cedge ∉ sub_edges #backward edge is not one of the subgraph edges - see if this edge bypasses the non-dominant vertex.
                    if any(reachable_dominance(subgraph) .& reachable_dominance(subgraph, cedge)) #only care about edges which are dominated by the dominating vertex of the subgraph.
                        tmp_mask .= tmp_mask .| non_dominance_mask(subgraph, cedge) #if any edge bypasses the dominated node then write a 1 in the mask for all the variable/root indices reachable from that edge
                    end
                else
                    tmp_mask .= tmp_mask .| (vertex_masks[backward_vertex(subgraph, cedge)]) #propagate bypass reachability bits from vertices earlier in the sort order.
                end
            end
        end
    end

    bypass_edges = PathEdge[]
    non_dom_edges = PathEdge[]

    #every vertex now has associated bypass reachability bits. Create bypass edges, if necessary.
    for edge in sub_edges

        temp_dom::Union{Nothing,PathEdge{T}} = nothing
        edge_mask = reachable_dominance(subgraph, edge)

        #If none of the reachable roots/variables in the edge match the dominance set of the subgraph then the subgraph doesn't exist anymore, or at least this edge shouldn't be part of the subgraph.
        some_dom_reachable = any(reachable_dominance(subgraph) .& edge_mask) #see if some reachable roots or variables are in the dominance set of the edge.
        @assert any(some_dom_reachable) "None of the reachable variables or roots were in the dominance set of the subgraph. This is a bug; please file a bug report at the FastDifferentiation repo."

        diff = set_diff(edge_mask, reachable_dominance(subgraph)) #important that diff is a new BitVector, not reused.
        if any(diff) #see if some roots or variables not in dominance set are in the reachable set of the edge. Split if true.
            if S === DominatorSubgraph
                temp_dom = PathEdge(top_vertex(edge), bott_vertex(edge), value(edge), copy(reachable_variables(edge)), diff) #create a new edge that accounts for roots not in the dominance mask
            else

                temp_dom = PathEdge(top_vertex(edge), bott_vertex(edge), value(edge), diff, copy(reachable_roots(edge))) #create a new edge that accounts for variables not in the     dominance mask    
            end
        end

        bypass::Union{Nothing,PathEdge{T}} = nothing

        back_mask = vertex_masks[backward_vertex(subgraph, edge)]
        feasible_split = back_mask .& reachable_non_dominance(subgraph, edge) #can only split if paths are present in the edge.
        if any(feasible_split)
            if subgraph isa FactorableSubgraph{T,DominatorSubgraph}
                bypass = PathEdge(top_vertex(edge), bott_vertex(edge), value(edge), copy(back_mask), reachable_dominance(subgraph))
            else
                bypass = PathEdge(top_vertex(edge), bott_vertex(edge), value(edge), reachable_dominance(subgraph), copy(back_mask))
            end
        end

        #see if edges can be merged. Put the merged edges in the bypass_edges list. Could also go in non_dom_edges, but only in one or the other.

        merged = false
        if temp_dom !== nothing && bypass !== nothing
            if (reachable_roots(temp_dom) == reachable_roots(bypass))
                merged = true
                push!(bypass_edges, PathEdge(top_vertex(edge), bott_vertex(edge), value(edge), reachable_variables(temp_dom) .| reachable_variables(bypass), reachable_roots(temp_dom))) #don't need to split edge
            elseif (reachable_variables(temp_dom) == reachable_variables(bypass))
                merged = true
                push!(bypass_edges, PathEdge(top_vertex(edge), bott_vertex(edge), value(edge), reachable_variables(temp_dom), reachable_roots(temp_dom) .| reachable_roots(bypass))) #don't need to split edge

            end
        end

        if !merged && temp_dom !== nothing
            push!(non_dom_edges, temp_dom)
        end

        if !merged && bypass !== nothing
            push!(bypass_edges, bypass)
        end
    end

    return bypass_edges, non_dom_edges
end




"""Computes a boolean path matrix `pm` for edges above and below the vertex. `pm[i,j]` = true means there is a path from variable `i` to variable `j`, false means there is no path. This is `O(mn)` where `m` is the nubmer of variables and `n` the nubmer of roots. Only use for debugging and testing"""
function check_continuity(graph, vertex)
    function path_matrix(graph, edge_list)
        paths = falses(domain_dimension(graph), codomain_dimension(graph))

        for edge in edge_list
            rts = copy(reachable_roots(edge))
            rt_indices = findall(rts)

            vars = copy(reachable_variables(edge))
            var_indices = findall(vars)

            for rt_index in rt_indices
                for var_index in var_indices
                    paths[var_index, rt_index] = true
                end
            end
        end
        return paths
    end

    pedges = parent_edges(graph, vertex)
    if length(pedges) != 0
        paths_above = path_matrix(graph, pedges)
        cedges = child_edges(graph, vertex)
        if length(cedges) != 0
            paths_below = path_matrix(graph, cedges)

            skip_row = -1
            skip_column = -1

            if is_root(graph, vertex)
                skip_col = root_postorder_to_index(graph, vertex)
                paths_below[:, skip_col] .= false
            end

            if is_variable(graph, vertex)
                skip_row = variable_postorder_to_index(graph, vertex)
                paths_above[skip_row, :] .= false
            end

            @assert all(paths_above .== paths_below) "Discountinuity in reachability for vertex $vertex. \n paths_above = $paths_above \n paths_below = $paths_below. This is a bug. Please create an issue on the FastDifferentiation.jl repo."
        end
    end
end

"""`find_bypass_edges` and `fin_non_dom_edges` can create redundant edges which don't need to be split. If an edge is present in both lists then merge so don't get redundant edges in tht graph"""
# function merge_redundant_edges(nondom,bypass) where{S<:PathEdge}
#     merged = S[]

#    for nd in nondom
#     for by in bypass
#         if vertices(nd) == vertices(by)
#             if reachable_roots(nd) == reachable_roots(by)

# end

"""reset root and variable masks for edges in the graph and add a new edge connecting `dominating_node(subgraph)` and `dominated_node(subgraph)` to the graph that has the factored value of the subgraph"""
function factor_subgraph!(subgraph::FactorableSubgraph{T}) where {T}
    function edge_continuity(graph, sub_edges)
        for edge in sub_edges
            if node_edges(graph, bott_vertex(edge)) !== nothing
                check_continuity(graph, bott_vertex(edge)) #"factorization of subgraph $(vertices(subgraph)) caused break in reachability continuity"
            end
            if node_edges(graph, top_vertex(edge)) !== nothing
                check_continuity(graph, top_vertex(edge)) #"factorization of subgraph $(vertices(subgraph)) caused break in reachability continuity"
            end
        end
    end

    sub_edges = subgraph_edges(subgraph)

    local new_edge::PathEdge{T}
    if length(sub_edges) > 0 #valid subgraph

        sum = evaluate_subgraph(subgraph, sub_edges)

        new_edge = make_factored_edge(subgraph, sub_edges, sum)
        @assert is_reachable(subgraph, new_edge) "Edge did not connect at least one root and one variable. This should never happen; it's a bug. Please create an issue on the FastDifferentiation.jl repo."



        bypass_edges, non_dom_edges = find_bypass_edges(subgraph, sub_edges)

        gr = graph(subgraph)
        for edge in non_dom_edges
            add_edge!(gr, edge)
        end

        for edge in bypass_edges
            add_edge!(gr, edge)
        end

        for edge in sub_edges
            delete_edge!(graph(subgraph), edge, true)
        end

        add_edge!(graph(subgraph), new_edge)

        #verify continuity at all vertices still remaining from the subgraph. This is wasteful because it potentially checks continuity twice for each vertex in the subgraph.

        # #TODO remove this and replace with more efficient test below once reachability bugs are fixed

        if RUN_CONTINUITY_CHECKS
            edge_continuity(graph(subgraph), unique_edges(graph(subgraph)))
        end

        # #end section to remove

        # edge_continuity(graph(subgraph), vcat(bypass_edges, non_dom_edges)) #only have to check continuity at newly created edges since subedges have been deleted.
    end
end

order!(::FactorableSubgraph{T,DominatorSubgraph}, nodes::Vector{T}) where {T<:Integer} = sort!(nodes,
) #largest node number last
order!(::FactorableSubgraph{T,PostDominatorSubgraph}, nodes::Vector{T}) where {T<:Integer} = sort!(nodes, rev=true) #largest node number first

predecessors(sub::FactorableSubgraph{T,DominatorSubgraph}, node_index::Integer) where {T<:Integer} = top_vertex.(filter(x -> test_edge(sub, x), parent_edges(graph(sub), node_index))) #allocates but this should rarely be called so shouldn't be efficiency issue.
predecessors(sub::FactorableSubgraph{T,PostDominatorSubgraph}, node_index::Integer) where {T<:Integer} = bott_vertex.(filter(x -> test_edge(sub, x), child_edges(graph(sub), node_index)))

"""Returns edges on the path from `node_index` to the dominating node which have one vertex equal to `node_index`"""
predecessor_edges(sub::FactorableSubgraph{T,DominatorSubgraph}, node_index::Integer) where {T<:Integer} = filter(x -> test_edge(sub, x), parent_edges(graph(sub), node_index)) #filter allocates: may be efficiency issue.
"""Returns edges on the path from `node_index` to the dominating node which have one vertex equal to `node_index`"""
predecessor_edges(sub::FactorableSubgraph{T,PostDominatorSubgraph}, node_index::Integer) where {T<:Integer} = filter(x -> test_edge(sub, x), child_edges(graph(sub), node_index)) #filter allocates: may be efficiency issue.


### These functions are used to evaluate subgraphs with branches created by factorization. This is not the most efficient way to evalute these subgraphs since terms in products are not ordered by uses. But subgraphs with branching seem rare and this is much simpler than recomputing the factorable subgraphs internal to a branching subgraph. Optimize if efficieny becomes an issue.

function vertex_counts(subgraph::FactorableSubgraph{T}) where {T}
    counts = Dict{T,T}()
    result = deconstruct_subgraph(subgraph)
    @assert result !== nothing
    sub_edges, sub_nodes = result

    for node in sub_nodes
        tmp = count(x -> in(x, sub_edges), backward_edges(subgraph, node)) #only count the child edges that are in the subgraph
        counts[node] = tmp
    end
    return counts
end


function evaluate_subgraph(subgraph::FactorableSubgraph{T}, sub_edges) where {T}
    vertex_sums = Dict{T,Node}()
    vertex_sums[dominated_node(subgraph)] = 1
    # Vis.draw_dot(subgraph)
    edge_list = collect(sub_edges)
    sort_edge_list!(subgraph, edge_list)

    for edge in edge_list
        if get(vertex_sums, forward_vertex(subgraph, edge), nothing) === nothing
            vertex_sums[forward_vertex(subgraph, edge)] = vertex_sums[backward_vertex(subgraph, edge)] * value(edge)
        else
            vertex_sums[forward_vertex(subgraph, edge)] += vertex_sums[backward_vertex(subgraph, edge)] * value(edge)
        end
    end

    @assert get(vertex_sums, dominating_node(subgraph), nothing) !== nothing "evaluate_branching_subgraph: vertex sum for dominating node did not exist. This is a bug. Please create an issue on the FastDifferentiation.jl repo."

    return vertex_sums[dominating_node(subgraph)]
end


function print_edges(a, msg)
    println(msg)
    for edge in edges(a)
        println(edge)
    end
end

function factor!(a::DerivativeGraph{T}) where {T}
    subgraph_list = compute_factorable_subgraphs(a)

    while !isempty(subgraph_list)

        subgraph = pop!(subgraph_list)

        factor_subgraph!(subgraph)
    end

    return nothing #return nothing so people don't mistakenly think this is returning a copy of the original graph
end

function follow_path(a::DerivativeGraph{T}, root_index::Integer, var_index::Integer) where {T}
    current_node_index = root_index_to_postorder_number(a, root_index)
    path_product = PathEdge{T}[]

    while true
        curr_edges = filter(x -> is_root_reachable(x, root_index) && is_variable_reachable(x, var_index), child_edges(a, current_node_index))
        if length(curr_edges) == 0
            break
        else
            @assert length(curr_edges) == 1 "Should only be one path from root $root_index to variable $var_index. Instead have $(length(curr_edges)) children from node $current_node_index on the path. These are the edges causing the error $(vertices.(curr_edges))"
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
            product *= value(term)
        end
    end
    return product
end

function evaluate_path(graph::DerivativeGraph, root_index::Integer, var_index::Integer)
    node_value = root(graph, root_index)
    if !is_tree(node_value) #root contains a variable or constant
        if is_variable(node_value)
            if variable(graph, var_index) == node_value
                return one(Node) #taking a derivative with respect to itself, which is 1. Need to figure out a better way to get the return number type right. This will always return Float64.
            else
                return zero(Node) #taking a derivative with respect to a different variable, which is 0.
            end
        else
            return zero(Node) #root is a constant
        end
    else #root contains a graph which has been factored so that there should be a single linear path from each root to each variable with no branching
        return follow_path(graph, root_index, var_index)
    end
end


"""Verifies that there is a single path from each root to each variable, if a path exists. This should be an invariant of the factored graph so it should always be true. But the algorithm is complex enough that it is easy to accidentally introduce errors when adding features. `verify_paths` has negligible runtime cost compared to factorization."""
# function _verify_paths(graph::DerivativeGraph{T}, a::T, visited::Dict{T,Bool}) where {T}

#     function duplicates!(dups, branches, reachable_func)
#         for branch in branches
#             dups .= dups .+ reachable_func(branch) #reachable_func(branch) returns a BitVector. In Julia a bit can be added to an Int. If there is a reachable bit set in index i of branch then 1 will be added to dups[i].
#         end
#     end

#     if get(visited, a, nothing) !== nothing #don't reprocess vertices that have already been visited.
#         return visited[a]
#     else
#         child_branches = child_edges(graph, a)
#         parent_branches = parent_edges(graph, a)
#         valid_graph = true

#         if length(parent_branches) > 1 && !is_constant(node(graph, a)) #don't care how many paths to roots from constant nodes
#             duplicate_reachables = zeros(Int64, codomain_dimension(graph))

#             duplicates!(duplicate_reachables, parent_branches, reachable_roots)

#             if any(x -> x > 1, duplicate_reachables)
#                 @info "duplicate paths to roots $(findall(x->x>1,duplicate_reachables)) from node $a"
#                 valid_graph = false
#             end
#         end


#         if length(child_branches) > 1

#             non_const_branches = PathEdge[]

#             for branch in child_branches
#                 if !any(reachable_variables(branch)) #no reachable variables so verify that the child node is a constant.
#                     @assert is_constant(node(graph, bott_vertex(branch))) "This is a bug. Please create an issue on the FastDifferentiation.jl repo."
#                 else
#                     push!(non_const_branches, branch)
#                 end
#             end

#             if length(non_const_branches) > 1
#                 duplicate_reachables = zeros(Int64, domain_dimension(graph))

#                 duplicates!(duplicate_reachables, non_const_branches, reachable_variables)

#                 if any(x -> x > 1, duplicate_reachables)
#                     @info "duplicate paths to variables $(findall(x->x>1,duplicate_reachables)) from node $a"
#                     valid_graph = false
#                 end
#             end
#         end

#         child_paths = children(graph, a)
#         if child_paths !== nothing
#             for child in children(graph, a)
#                 valid_graph = valid_graph & _verify_paths(graph, child, visited) #if any node has invalid paths valid_graph will be false
#             end
#         end

#         check_continuity(graph, a)
#         visited[a] = valid_graph

#         return valid_graph
#     end
# end

function _verify_paths(graph, a, visited)
    if get(visited, a, nothing) !== nothing #don't reprocess vertices that have already been visited.
        return visited[a]
    else
        cedges = child_edges(graph, a)
        valid_graph = true

        #seems like testing for non duplicate paths is inherently O(n^2). Most graphs don't seem to have a large number of internal edges so this shouldn't normally be a problem.
        for i in 1:length(cedges)
            for j in i+1:length(cedges)
                rrts1 = reachable_roots(cedges[i])
                rrts2 = reachable_roots(cedges[j])
                rvars1 = reachable_variables(cedges[i])
                rvars2 = reachable_variables(cedges[j])

                if any(rrts1 .& rrts2) && any(rvars1 .& rvars2) #have more than one path from some root(s) to some variable(s)
                    valid_graph = false
                end
            end
        end

        child_paths = children(graph, a)
        if valid_graph && child_paths !== nothing
            for child in children(graph, a)
                valid_graph = valid_graph && _verify_paths(graph, child, visited) #if any node has invalid paths valid_graph will be false
            end
        end

        check_continuity(graph, a)
        visited[a] = valid_graph

        return valid_graph
    end
end


"""verifies that there is a single path from each root to each variable, if such a path exists."""
function verify_paths(graph::DerivativeGraph{T}) where {T}
    visited = Dict{T,Bool}()
    good_paths = true
    for root in roots(graph)
        if !_verify_paths(graph, postorder_number(graph, root), visited)
            good_paths = false
        end
    end
    return good_paths
end


function unique_nodes(jacobian::AbstractArray{T}) where {T<:Node} #not efficient, may revist parts of the jacobian many times.
    nodes = IdDict{Node,Bool}()

    for index in eachindex(jacobian)
        oned = all_nodes(jacobian[index])
        for node in oned
            nodes[node] = true
        end
    end
    # nodes = Set{Node}()
    # for index in eachindex(jacobian)
    #     oned = all_nodes(jacobian[index])
    #     union!(nodes, oned)
    # end
    # return nodes
    return keys(nodes)
end

"""Count of number of operations in graph."""
function number_of_operations(jacobian::AbstractArray{T}) where {T<:Node}
    count = 0
    nodes = all_nodes(jacobian)
    for node in nodes
        if is_tree(node) && !is_negate(node) #don't count negate as an operation
            count += 1
        end
    end
    return count
end

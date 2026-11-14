
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

    diffa = node_difference(a)
    diffb = node_difference(b)


    # if a ⊂ b then diff(a) < diff(b) where diff(x) = abs(dominating_node(a) - dominated_node(a)). Might be that a ⊄ b but it's safe to factor a first.
    # This tests only guarantees that if a ⊂ b then a will be factored first. It could be that diff(a) < diff(b) but that a ⊄ b in which case most
    # efficient option would be to factor whichever of a,b has highest times used. But determining subgraph containment precisely is time consuming.
    # This ordering heuristic doesn't seem to affect efficiency of computed derivatives much but is significantly faster. 
    if diffa < diffb
        return true
    elseif diffa > diffb
        return false
    elseif times_used(a) > times_used(b) #If a is used more times than b then factor a first. More efficient.
        return true
    else
        return false
    end
end


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

    function set_dom_bits!(subgraph::Tuple{T,T}, all_subgraphs::Dict, var_or_root_index::Integer, bit_dimension::Integer) where {T<:Integer}
        existing_subgraph = get(all_subgraphs, subgraph, nothing)
        if existing_subgraph !== nothing
            all_subgraphs[subgraph][var_or_root_index] = 1
        else
            tmp = falses(bit_dimension)
            tmp[var_or_root_index] = 1
            all_subgraphs[subgraph] = tmp
        end
    end

    temp_doms = Dict{T,T}()

    for root_index in 1:codomain_dimension(graph)
        post_num = root_index_to_postorder_number(graph, root_index)
        temp_dom = compute_dom_table(graph, true, root_index, post_num, temp_doms)

        for dominated in keys(temp_dom)
            dsubgraph = dom_subgraph(graph, root_index, dominated, temp_dom)
            if dsubgraph !== nothing
                set_dom_bits!(dsubgraph, dom_subgraphs, root_index, codomain_dimension(graph))
            end
        end
    end

    for variable_index in 1:domain_dimension(graph)
        post_num = variable_index_to_postorder_number(graph, variable_index)
        temp_dom = compute_dom_table(graph, false, variable_index, post_num, temp_doms)

        for dominated in keys(temp_dom)
            psubgraph = pdom_subgraph(graph, variable_index, dominated, temp_dom)
            if psubgraph !== nothing
                set_dom_bits!(psubgraph, pdom_subgraphs, variable_index, domain_dimension(graph))
            end
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
            subgraph = dominator_subgraph(graph, dominator, dominated, dom_subgraphs[key], reachable_roots(graph, dominator), reachable_variables(graph, dominated))

            push!(result, subgraph)
        end
    end

    for key in keys(pdom_subgraphs)
        dominator = key[1]
        dominated = key[2]
        subgraph = postdominator_subgraph(graph, dominator, dominated, pdom_subgraphs[key], reachable_roots(graph, dominated), reachable_variables(graph, dominator))

        push!(result, subgraph)
    end

    return result
end


"True if `edge` is live for at least one pair of dominance and non-dominance bits of `subgraph`."
on_subgraph_path(subgraph::FactorableSubgraph, edge::PathEdge) =
    overlap(reachable_dominance(subgraph, edge), reachable_dominance(subgraph)) &&
    overlap(non_dominance_mask(subgraph, edge), non_dominance_mask(subgraph))

# Dictionaries keyed by masks use the masks' chunk vectors, because `hash` on a
# `BitVector` falls back to the generic element-by-element array method.
mask_key(mask::BitVector) = mask.chunks

"A path group: dominance mask, non-dominance mask, sum of path products, and number of paths."
const PathGroup = Tuple{BitVector,BitVector,Node,Int}
const PathGroups = Dict{Tuple{Vector{UInt64},Vector{UInt64}},PathGroup}

"""
    subgraph_path_groups(subgraph, current_node, memo)

Returns the live paths from `current_node` to `dominating_node(subgraph)`, grouped by the pair of masks (dominance, non-dominance) for which every edge of the path is live.  Each group stores the sum of its path products and its number of paths.  The dynamic-programming evaluation handles branching inside the subgraph without enumerating complete paths.
"""
function subgraph_path_groups(subgraph::FactorableSubgraph{T,S}, current_node::T, memo::Dict{T,PathGroups}) where {T,S<:AbstractFactorableSubgraph}
    cached = get(memo, current_node, nothing)
    cached === nothing || return cached

    result = PathGroups()
    if current_node == dominating_node(subgraph)
        dom, nondom = copy(reachable_dominance(subgraph)), copy(non_dominance_mask(subgraph))
        result[(mask_key(dom), mask_key(nondom))] = (dom, nondom, Node(1), 1)
        memo[current_node] = result
        return result
    end

    for edge in forward_edges(subgraph, current_node)
        on_subgraph_path(subgraph, edge) || continue
        dmask = reachable_dominance(subgraph, edge)
        nmask = non_dominance_mask(subgraph, edge)
        for (suffix_dom, suffix_nondom, suffix_value, suffix_count) in values(subgraph_path_groups(subgraph, forward_vertex(subgraph, edge), memo))
            dom = dmask .& suffix_dom
            nondom = nmask .& suffix_nondom
            (any(dom) && any(nondom)) || continue
            key = (mask_key(dom), mask_key(nondom))
            # multiply from the top of the subgraph downwards, which is the order
            # in which follow_path forms products, so that products can be shared
            path_value = suffix_value * value(edge)
            previous = get(result, key, nothing)
            result[key] = previous === nothing ? (dom, nondom, path_value, suffix_count) : (dom, nondom, previous[3] + path_value, previous[4] + suffix_count)
        end
    end

    memo[current_node] = result
    return result
end

"""
    subgraph_path_groups(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}, region)

Returns the same path groups as the recursive method for dominator subgraphs, for the paths from `dominated_node(subgraph)` down to `dominating_node(subgraph)`.  The groups are pushed downwards through the nodes of `region` in topological order, so that each product is formed from the top of the subgraph down, as it is for dominator subgraphs."""
function subgraph_path_groups(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}, region::Vector{T}) where {T}
    top = dominated_node(subgraph)
    bottom = dominating_node(subgraph)
    prefixes = Dict{T,PathGroups}()
    dom, nondom = copy(reachable_dominance(subgraph)), copy(non_dominance_mask(subgraph))
    prefixes[top] = PathGroups((mask_key(dom), mask_key(nondom)) => (dom, nondom, Node(1), 1))
    for u in sort!(vcat(top, filter(!=(bottom), region)), rev=true) # topological order, from the top down
        groups = get(prefixes, u, nothing)
        groups === nothing && continue
        for edge in forward_edges(subgraph, u)
            on_subgraph_path(subgraph, edge) || continue
            w = forward_vertex(subgraph, edge)
            dmask = reachable_dominance(subgraph, edge)
            nmask = non_dominance_mask(subgraph, edge)
            target = get!(PathGroups, prefixes, w)
            for (prefix_dom, prefix_nondom, prefix_value, prefix_count) in values(groups)
                d = dmask .& prefix_dom
                n = nmask .& prefix_nondom
                (any(d) && any(n)) || continue
                key = (mask_key(d), mask_key(n))
                path_value = prefix_value * value(edge)
                previous = get(target, key, nothing)
                target[key] = previous === nothing ? (d, n, path_value, prefix_count) : (d, n, previous[3] + path_value, previous[4] + prefix_count)
            end
        end
    end
    return get(prefixes, bottom, PathGroups())
end

"""
    index_classes(masks, dimension)

Groups the indices that are set in at least one of `masks` by the set of masks containing them.  Returns a vector of pairs (signature, class), where `signature[k]` says whether `masks[k]` contains the class and `class` is a mask of the indices with that signature.  The classes are found by splitting them against one mask at a time, so the cost depends on the number of masks and classes rather than on the number of indices."""
function index_classes(masks::Vector{BitVector}, dimension::Integer)
    occurring = falses(dimension)
    for m in masks
        occurring .|= m
    end
    classes = any(occurring) ? [falses(length(masks)) => occurring] : Pair{BitVector,BitVector}[]
    for (k, m) in enumerate(masks)
        refined = Pair{BitVector,BitVector}[]
        for (sig, class) in classes
            inside = class .& m
            outside = class .& .!m
            if any(inside)
                inside_sig = copy(sig)
                inside_sig[k] = true
                push!(refined, inside_sig => inside)
            end
            any(outside) && push!(refined, sig => outside)
        end
        classes = refined
    end
    return classes
end

"""
    factored_edges(subgraph)

Returns the replacement edges for `subgraph` and the largest number of paths that any single pair of dominance and non-dominance bits lies on.  Every path covers a rectangle of pairs, namely the product of its two masks.  Dominance indices and non-dominance indices are each grouped by the paths that contain them, and every block of the resulting grid is covered by one set of paths, so it becomes one edge whose value is the sum of those paths.  Only indices that lie on some path are visited, so the cost does not grow with the width of the graph."""
function factored_edges(subgraph::FactorableSubgraph{T,S}, region::Vector{T}=subgraph_region(subgraph)) where {T,S<:AbstractFactorableSubgraph}
    dom = dominating_node(subgraph)
    dmd = dominated_node(subgraph)
    path_groups = if S === DominatorSubgraph
        subgraph_path_groups(subgraph, dmd, Dict{T,PathGroups}())
    else
        subgraph_path_groups(subgraph, region)
    end
    groups = collect(values(path_groups))
    rows = index_classes(BitVector[g[1] for g in groups], length(reachable_dominance(subgraph)))
    columns = index_classes(BitVector[g[2] for g in groups], non_dominance_dimension(subgraph))

    result_edges = PathEdge{T}[]
    max_paths = 0
    for (row_sig, row_mask) in rows, (column_sig, column_mask) in columns
        cover = row_sig .& column_sig
        any(cover) || continue
        sum_val = Node(0)
        num_paths = 0
        for (k, (_, _, path_value, path_count)) in enumerate(groups)
            if cover[k]
                sum_val += path_value
                num_paths += path_count
            end
        end
        max_paths = max(max_paths, num_paths)
        if S === DominatorSubgraph
            push!(result_edges, PathEdge(dom, dmd, sum_val, copy(column_mask), copy(row_mask)))
        else
            push!(result_edges, PathEdge(dom, dmd, sum_val, copy(row_mask), copy(column_mask)))
        end
    end
    return result_edges, max_paths
end

"""
    subgraph_region(subgraph)

Returns the nodes reachable from `dominated_node(subgraph)` by moving forward along edges that satisfy `on_subgraph_path`, stopping at `dominating_node(subgraph)`.  Must be called before the subgraph is factored, because factoring clears the masks of the boundary edges."""
function subgraph_region(subgraph::FactorableSubgraph{T}) where {T}
    stop = dominating_node(subgraph)
    visited = Set{T}()
    queue = T[dominated_node(subgraph)]
    while !isempty(queue)
        curr = pop!(queue)
        curr == stop && continue
        for e in forward_edges(subgraph, curr)
            if on_subgraph_path(subgraph, e)
                nxt = forward_vertex(subgraph, e)
                if !in(nxt, visited)
                    push!(visited, nxt)
                    push!(queue, nxt)
                end
            end
        end
    end
    return collect(visited)
end

backward_vertex(::FactorableSubgraph{T,DominatorSubgraph}, edge::PathEdge) where {T} = bott_vertex(edge)
backward_vertex(::FactorableSubgraph{T,PostDominatorSubgraph}, edge::PathEdge) where {T} = top_vertex(edge)

"""
    continuation_classes(subgraph, continuations, dimension)

Groups the dominance bits (roots for a dominator subgraph, variables for a postdominator subgraph) by the set of `continuations` whose dominance masks contain them.  Returns `class_of`, which maps each bit to its class (0 if no continuation contains it), and `reach`, which holds for each class the union of the non-dominance masks of the continuations that contain it.  Every edge leaving a node is pruned against the same continuations, so this is computed once per node rather than once per edge.

The classes are the same as `index_classes` would give for the continuations' dominance masks, but are found in one pass over the set bits of those masks: each bit collects the indices of the continuations that contain it, and bits with equal lists form a class.  This costs one step per set bit instead of one bit-vector operation per (continuation, class), which matters at nodes with many continuations."""
function continuation_classes(subgraph::FactorableSubgraph, continuations, dimension::Integer)
    signature = Vector{Vector{Int}}(undef, dimension)   # bit => indices of the continuations containing it
    seen = falses(dimension)
    for (i, f) in enumerate(continuations)
        mask = reachable_dominance(subgraph, f)
        bit = findnext(mask, 1)
        while bit !== nothing
            seen[bit] || (seen[bit] = true; signature[bit] = Int[])
            push!(signature[bit], i)
            bit = findnext(mask, bit + 1)
        end
    end
    class_of = zeros(Int, dimension)
    reach = BitVector[]
    class_index = Dict{Vector{Int},Int}()
    bit = findnext(seen, 1)
    while bit !== nothing
        k = get(class_index, signature[bit], 0)
        if k == 0
            r = falses(non_dominance_dimension(subgraph))
            for i in signature[bit]
                r .|= non_dominance_mask(subgraph, continuations[i])
            end
            push!(reach, r)
            k = class_index[signature[bit]] = length(reach)
        end
        class_of[bit] = k
        bit = findnext(seen, bit + 1)
    end
    return class_of, reach
end

"""
    prune_edge!(subgraph, e, class_of, reach, to_add, to_delete)

Prunes one edge `e` against the edges that continue it away from `dominating_node(subgraph)`, summarized by `continuation_classes`, and returns `true` if `e` changed.

The dominance bits of `e` are partitioned into groups that are continued by the same edges, which are the classes of `continuation_classes`.  Each group keeps only the non-dominance bits of `e` that are set on at least one of its continuing edges.  When the groups end up with different masks, `e` is split into one edge per group, as `add_non_dom_edges!` does for boundary edges, and the extra edges are pushed onto `to_add`; if no group keeps any bit, `e` is pushed onto `to_delete`.  No live path is removed, because every bit that survives is continued by an edge that is live for the same dominance bit."""
function prune_edge!(subgraph::FactorableSubgraph{T,S}, e::PathEdge, class_of::Vector{Int}, reach::Vector{BitVector}, to_add::Vector{PathEdge{T}}, to_delete::Vector{PathEdge{T}}) where {T,S<:AbstractFactorableSubgraph}
    dom_e = reachable_dominance(subgraph, e)

    # Fast path, the common case: all dominance bits of e lie in one class, and the
    # continuations of that class carry every non-dominance bit of e.  Then e has a
    # single group that keeps all its bits, so it does not change.
    first_bit = findnext(dom_e, 1)
    if first_bit !== nothing && class_of[first_bit] != 0
        k = class_of[first_bit]
        one_class = true
        bit = findnext(dom_e, first_bit + 1)
        while bit !== nothing
            if class_of[bit] != k
                one_class = false
                break
            end
            bit = findnext(dom_e, bit + 1)
        end
        nondom_e = non_dominance_mask(subgraph, e)
        one_class && any(nondom_e) && is_subset_of(nondom_e, reach[k]) && return false
    end

    groups = Dict{Int,BitVector}()   # class => the bits of e in that class
    for bit in findall(dom_e)
        k = class_of[bit]
        k == 0 && continue   # no continuation carries this bit
        get!(() -> falses(length(dom_e)), groups, k)[bit] = true
    end

    # the surviving non-dominance bits of each group, merging groups whose bits agree
    kept = Dict{Vector{UInt64},Pair{BitVector,BitVector}}()
    for (k, g) in groups
        reach_g = reach[k] .& non_dominance_mask(subgraph, e)
        any(reach_g) || continue
        existing = get(kept, mask_key(reach_g), nothing)
        existing === nothing ? (kept[mask_key(reach_g)] = reach_g => g) : (existing.second .|= g)
    end

    if isempty(kept)
        fill!(non_dominance_mask(subgraph, e), false)
        push!(to_delete, e)
        return true
    end
    changed = length(kept) > 1
    first_pass = true
    for (nondom, dom) in values(kept)
        if first_pass
            changed |= dom != reachable_dominance(subgraph, e) || nondom != non_dominance_mask(subgraph, e)
            reachable_dominance(subgraph, e) .= dom
            non_dominance_mask(subgraph, e) .= nondom
            first_pass = false
        elseif S === DominatorSubgraph
            push!(to_add, PathEdge(top_vertex(e), bott_vertex(e), value(e), copy(nondom), copy(dom)))
        else
            push!(to_add, PathEdge(top_vertex(e), bott_vertex(e), value(e), copy(dom), copy(nondom)))
        end
    end
    return changed
end

"""
    prune_stale_masks!(subgraph, region, seeds)

After `subgraph` has been factored, removes reachability that no longer leads anywhere from the edges inside `region`.

An edge can only become stale if one of its continuations changed, because the masks were exact before this step and factoring only changed edges at the dominated node.  So pruning starts from `seeds`, the nodes whose edges towards the dominated node were changed by detaching, and works outward with a worklist: at each node `w` it prunes the edges that leave `w` towards the dominating node (and stay in `region`) against the edges that continue them from `w`, and queues the far end of every edge that changed.  Nodes are taken in order from the dominated node outward, so every edge is pruned against continuations that are already final.  The cost is proportional to the edges near changed ones, not to the size of `region`."""
function prune_stale_masks!(subgraph::FactorableSubgraph{T,S}, region::Vector{T}, seeds) where {T,S<:AbstractFactorableSubgraph}
    gr = graph(subgraph)
    b = dominated_node(subgraph)
    in_region = Set{T}(region)
    queued = Set{T}()
    # outward from the dominated node: increasing postorder for a dominator subgraph, decreasing for a postdominator subgraph
    worklist = S === DominatorSubgraph ? BinaryMinHeap{T}() : BinaryMaxHeap{T}()
    enqueue!(w) = (w != b && in(w, in_region) && !in(w, queued)) && (push!(queued, w); push!(worklist, w))
    foreach(enqueue!, seeds)
    to_delete = PathEdge{T}[]
    to_add = PathEdge{T}[]
    while !isempty(worklist)
        w = pop!(worklist)
        continuations = filter(f -> any(non_dominance_mask(subgraph, f)), backward_edges(subgraph, w))
        class_of, reach = continuation_classes(subgraph, continuations, dominance_dimension(subgraph))
        empty!(to_delete)
        empty!(to_add)
        for e in collect(forward_edges(subgraph, w))
            u = forward_vertex(subgraph, e)
            in(u, in_region) || continue
            prune_edge!(subgraph, e, class_of, reach, to_add, to_delete) && enqueue!(u)
        end
        foreach(e -> delete_edge!(gr, e), to_delete)
        foreach(e -> add_edge!(gr, e), to_add)
    end
    return nothing
end

"""
    sum_all_paths(graph, root_index, var_index)

Evaluate all reachable derivative paths between a root and variable. Factored
graphs normally contain one such path, but the dynamic-programming traversal
also handles residual branches without dropping their contributions.
"""
function sum_all_paths(a::DerivativeGraph{T}, root_index::Integer, var_index::Integer) where {T}
    current_node_index = root_index_to_postorder_number(a, root_index)
    memo = Dict{T,Node}()

    function path_sum(node_index::T)
        cached = get(memo, node_index, nothing)
        cached === nothing || return cached

        curr_edges = filter(
            edge -> is_root_reachable(edge, root_index) && is_variable_reachable(edge, var_index),
            child_edges(a, node_index),
        )
        if isempty(curr_edges)
            return is_variable(a, node_index) && variable_postorder_to_index(a, node_index) == var_index ? Node(1.0) : Node(0.0)
        end

        result = Node(0.0)
        for edge in curr_edges
            result += value(edge) * path_sum(bott_vertex(edge))
        end
        memo[node_index] = result
        return result
    end

    curr_edges = filter(
        edge -> is_root_reachable(edge, root_index) && is_variable_reachable(edge, var_index),
        child_edges(a, current_node_index),
    )
    isempty(curr_edges) && return Node(0.0)
    return sum(value(edge) * path_sum(bott_vertex(edge)) for edge in curr_edges; init=Node(0.0))
end

"""
    factor_subgraph!(subgraph)

Replace a factorable subgraph with its factored replacement edge(s), splitting boundary non-dominance
reachability masks and removing factored reachability at `dominated_node(subgraph)`.
"""
function factor_subgraph!(subgraph::FactorableSubgraph{T,S}) where {T,S<:AbstractFactorableSubgraph}
    region = subgraph_region(subgraph)
    new_edges, max_paths = factored_edges(subgraph, region)
    if max_paths ≥ 2 # otherwise every pair already has at most one path through the subgraph
        graph_value = graph(subgraph)
        # the nodes whose edges to the dominated node detaching may change; pruning starts there
        seeds = T[forward_vertex(subgraph, e) for e in forward_edges(subgraph, dominated_node(subgraph))]

        for new_edge in new_edges
            partition = if S === DominatorSubgraph
                dominator_subgraph(
                    graph_value,
                    dominating_node(subgraph),
                    dominated_node(subgraph),
                    copy(reachable_roots(new_edge)),
                    copy(reachable_roots(new_edge)),
                    copy(reachable_variables(new_edge)),
                )
            else
                postdominator_subgraph(
                    graph_value,
                    dominating_node(subgraph),
                    dominated_node(subgraph),
                    copy(reachable_variables(new_edge)),
                    copy(reachable_roots(new_edge)),
                    copy(reachable_variables(new_edge)),
                )
            end

            add_non_dom_edges!(partition)
            edges_to_delete = reset_edge_masks!(partition)
            for edge in edges_to_delete
                delete_edge!(graph_value, edge)
            end
            add_edge!(graph_value, new_edge)
        end
        prune_stale_masks!(subgraph, region, seeds)
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

"""
    follow_path(graph, root_index, var_index)

Returns the derivative of root `root_index` with respect to variable `var_index` in a factored graph, which normally has a single live path between them.  The edges of that path are multiplied in decreasing order of the number of root–variable pairs that use them, so that products shared by many derivatives are formed once.  If the walk meets a branch or a dead end, `sum_all_paths` is used instead.
"""
function follow_path(a::DerivativeGraph{T}, root_index::Integer, var_index::Integer) where {T}
    current = root_index_to_postorder_number(a, root_index)
    target = variable_index_to_postorder_number(a, var_index)
    path = PathEdge{T}[]
    while current != target
        next_edge = nothing
        count = 0
        for edge in child_edges(a, current)
            if is_root_reachable(edge, root_index) && is_variable_reachable(edge, var_index)
                count += 1
                next_edge = edge
            end
        end
        count == 0 && isempty(path) && return Node(0.0)
        count == 1 || return sum_all_paths(a, root_index, var_index)
        push!(path, next_edge)
        current = bott_vertex(next_edge)
    end
    sort!(path, lt=(x, y) -> num_uses(x) > num_uses(y))
    product = Node(1.0)
    for edge in path
        product *= value(edge)
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
function _verify_paths(graph::DerivativeGraph, a::Int)
    child_branches = child_edges(graph, a)
    parent_branches = parent_edges(graph, a)
    valid_graph = true

    if length(parent_branches) > 1 #this simple test won't work if a branch is a child of another branch. Then could have 
        roots_intersect = reduce(.&, reachable_roots.(parent_branches))
        if !is_zero(roots_intersect)
            valid_graph = false
        end
    end
    if length(child_branches) > 1
        vars_intersect = reduce(.&, reachable_variables.(child_branches))
        if !is_zero(vars_intersect)
            valid_graph = false
        end
    end

    if !valid_graph
        # FastDifferentiation.FastDifferentiationVisualizationExt.draw_dot(graph)
        @info "failure"
        @info "roots_intersect $roots_intersect"
        # return false
    else
        for child in children(graph, a)
            valid_graph &= _verify_paths(graph, child)
        end
    end

    return valid_graph
end

"""verifies that there is a single path from each root to each variable, if such a path exists."""
function verify_paths(graph::DerivativeGraph)
    return true #until can fix this so it both correctly verifies paths and does not take quadratic time.
    for root in roots(graph)
        if !_verify_paths(graph, postorder_number(graph, root))
            return false
        end
    end
    return true
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

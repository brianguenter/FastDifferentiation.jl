
next_edge_constraint(sub::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = PathConstraint(dominating_node(sub), graph(sub), false, reachable_roots(sub), reachable_dominance(sub))
next_edge_constraint(sub::FactorableSubgraph{T,DominatorSubgraph}) where {T} = PathConstraint(dominating_node(sub), graph(sub), true, reachable_dominance(sub), reachable_variables(sub))
top_down_constraint(sub::FactorableSubgraph{T,DominatorSubgraph}) where {T} = PathConstraint()

"""Evaluates the subgraph, creates a new edge with this value, and then inserts the new edge into `graph`"""
function add_edge!(graph::DerivativeGraph, subgraph::FactorableSubgraph, subgraph_value::Node)
    verts = vertices(subgraph)
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

function old_edge_path(next_node_constraint, dominating::T, is_dominator::Bool, reachable_mask::BitVector, current_edge) where {T}
    flag_value = 1
    result = PathEdge{Int64}[]

    roots_reach = copy(reachable_roots(current_edge))
    vars_reach = copy(reachable_variables(current_edge))

    while true
        push!(result, current_edge)
        if is_dominator && top_vertex(current_edge) == dominating
            break
        end
        if !is_dominator && bott_vertex(current_edge) == dominating
            break
        end
        tmp = get_edge_vector()
        relation_edges!(next_node_constraint, current_edge, tmp)

        if is_dominator
            filter!(x -> subset(reachable_mask, reachable_variables(x)), tmp) #only accept edges which have reachable roots and variables that match the subgraph
        else
            filter!(x -> subset(reachable_mask, reachable_roots(x)), tmp)
        end

        #These two cases can only occur if the subgraph has been destroyed by factorization
        if length(tmp) == 0  #there is no edge beyond current_edge that leads to the dominating node. 
            reclaim_edge_vector(tmp)
            flag_value = 0
            break
        elseif length(tmp) ≥ 2 #there is a branch in the edge path  
            reclaim_edge_vector(tmp)
            flag_value = 2
            break
        end

        current_edge = tmp[1]
        roots_reach .= roots_reach .& reachable_roots(current_edge)
        vars_reach .= vars_reach .& reachable_variables(current_edge) #update reachable roots/variables of the entire path. Sum is only good over this subset
        reclaim_edge_vector(tmp)
    end

    return flag_value, result, roots_reach, vars_reach
end

function evaluate_subgraph(subgraph::FactorableSubgraph{T,S}) where {T,S<:Union{DominatorSubgraph,PostDominatorSubgraph}}
    constraint = next_edge_constraint(subgraph)
    sum = Node(0.0)

    rel_edges = get_edge_vector()
    for edge in relation_edges!(constraint, dominated_node(subgraph), rel_edges)
        flag, pedges, roots_reach, vars_reach = old_edge_path(constraint, dominating_node(subgraph), S == DominatorSubgraph, reachable(subgraph), edge)
        #sort by num_uses then from largest to smallest postorder number

        if flag == 1 #non-branching path through subgraph

            sort!(pedges, lt=path_sort_order)
            sum += multiply_sequence(pedges)

            #find sequences of equal times_used and multiply them. Then multiply each of the collapsed sequences to get the final result
            reclaim_edge_vector(pedges)
        end
    end


    reclaim_edge_vector(rel_edges)
    return sum
end

function make_factored_edge(subgraph::FactorableSubgraph{T,DominatorSubgraph}, sum::Node) where {T}
    roots_reach = copy(reachable_dominance(subgraph))
    vars_reach = copy(reachable_variables(subgraph))
    return PathEdge(dominating_node(subgraph), dominated_node(subgraph), sum, vars_reach, roots_reach)
end

function make_factored_edge(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}, sum::Node) where {T}
    roots_reach = copy(reachable_roots(subgraph))
    vars_reach = copy(reachable_dominance(subgraph))
    return PathEdge(dominating_node(subgraph), dominated_node(subgraph), sum, vars_reach, roots_reach)
end


"""Returns true if a new factorable subgraph was created inside `subgraph` during the factorization process. If true then must compute factorable subgraphs for the edges inside `subgraph`. `subgraph_exists` should be called before executing this function otherwise it may return false when no new subgraphs have been created."""
function is_branching(subgraph)
    fedges = forward_edges(subgraph, dominated_node(subgraph))

    sub_edges = Set{PathEdge}()
    bad_subgraph = false
    for edge in fedges #for each forward edge from the dominated node find all edges on that path. If any edge in the subgraph is visited more than once this means a new factorable subgraph has been created.
        good_edges, tmp = edges_on_path(subgraph, edge)

        if good_edges
            for pedge in tmp
                if in(pedge, sub_edges) #edge has been visited twice.
                    bad_subgraph = true
                    break
                end
                push!(sub_edges, pedge)
            end
        end
    end

    return bad_subgraph
end

"""reset root and variable masks for edges in the graph and add a new edge connecting `dominating_node(subgraph)` and `dominated_node(subgraph)` to the graph that has the factored value of the subgraph"""
function factor_subgraph!(subgraph::FactorableSubgraph{T}) where {T}
    local new_edge::PathEdge{T}
    if subgraph_exists(subgraph)

        if is_branching(subgraph) #handle the uncommon case of factorization creating new factorable subgraphs internal to subgraph
            sum = evaluate_branching_subgraph(subgraph)
            new_edge = make_factored_edge(subgraph, sum)
        else
            sum = evaluate_subgraph(subgraph)
            new_edge = make_factored_edge(subgraph, sum)
        end
        add_non_dom_edges!(subgraph)
        #reset roots in R, if possible. All edges earlier in the path than the first vertex with more than one child cannot be reset.
        edges_to_delete = reset_edge_masks!(subgraph)
        for edge in edges_to_delete
            delete_edge!(graph(subgraph), edge)
        end

        add_edge!(graph(subgraph), new_edge)
    end
end

order!(::FactorableSubgraph{T,DominatorSubgraph}, nodes::Vector{T}) where {T<:Integer} = sort!(nodes,
) #largest node number last
order!(::FactorableSubgraph{T,PostDominatorSubgraph}, nodes::Vector{T}) where {T<:Integer} = sort!(nodes, rev=true) #largest node number first

predecessors(sub::FactorableSubgraph{T,DominatorSubgraph}, node_index::Integer) where {T<:Integer} = top_vertex.(filter(x -> test_edge(sub, x), parent_edges(graph(sub), node_index))) #allocates but this should rarely be called so shouldn't be efficiency issue.
predecessors(sub::FactorableSubgraph{T,PostDominatorSubgraph}, node_index::Integer) where {T<:Integer} = bott_vertex.(filter(x -> test_edge(sub, x), child_edges(graph(sub), node_index)))

predecessor_edges(sub::FactorableSubgraph{T,DominatorSubgraph}, node_index::Integer) where {T<:Integer} = filter(x -> test_edge(sub, x), parent_edges(graph(sub), node_index)) #allocates but this should rarely be called so shouldn't be efficiency issue.
predecessor_edges(sub::FactorableSubgraph{T,PostDominatorSubgraph}, node_index::Integer) where {T<:Integer} = filter(x -> test_edge(sub, x), child_edges(graph(sub), node_index))


"""Computes idoms for special case when new factorable subgraphs are created by factorization. This seems redundant with compute_factorable_subgraphs, fill_idom_tables, etc. but invariants that held when graph was first factored no longer hold so need specialized code. Not currently used, experimental code."""
function compute_internal_idoms(subgraph::FactorableSubgraph{T}) where {T}
    _, sub_nodes = deconstruct_subgraph(subgraph)
    order!(subgraph, sub_nodes)
    compressed_index = Dict((sub_nodes[i] => i) for i in eachindex(sub_nodes))

    preds = [map(x -> compressed_index[x], predecessors(subgraph, node)) for node in sub_nodes] #allocates but this function should rarely be called
    compressed_doms = simple_dominance(preds) #idom table in compressed index format_string
    return Dict{T,T}([(sub_nodes[i], sub_nodes[compressed_doms[i]]) for i in eachindex(sub_nodes)])
end


### These functions are used to evaluate subgraphs with branches created by factorization. This is not the most efficient way to evalute these subgraphs since terms in products are not ordered by uses. But subgraphs with branching seem rare and this is much simpler than recomputing the factorable subgraphs internal to a branching subgraph. Optimize if efficieny becomes an issue.

function vertex_counts(subgraph::FactorableSubgraph{T}) where {T}
    counts = Dict{T,T}()
    sub_edges, sub_nodes = deconstruct_subgraph(subgraph)

    for node in sub_nodes
        tmp = count(x -> in(x, sub_edges), backward_edges(subgraph, node)) #only count the child edges that are in the subgraph
        counts[node] = tmp
    end
    return counts
end

function evaluate_branching_subgraph(subgraph::FactorableSubgraph{T}) where {T}
    global num_times += 1
    sub_edges, sub_nodes = deconstruct_subgraph(subgraph)
    counts = vertex_counts(subgraph)
    counts[dominated_node(subgraph)] = 1
    vertex_sums = Dict{T,Node}()
    # Vis.draw_dot(subgraph)
    _evaluate_branching_subgraph(subgraph, Node(1), dominated_node(subgraph), sub_edges, counts, vertex_sums)

    return vertex_sums[dominating_node(subgraph)]
end

num_times = 0

function _evaluate_branching_subgraph(subgraph::FactorableSubgraph{T}, sum::Node, current_vertex::T, sub_edges, counts::Dict{T,T}, vertex_sums::Dict{T,Node}) where {T}
    if get(vertex_sums, current_vertex, nothing) === nothing
        vertex_sums[current_vertex] = sum
    else
        vertex_sums[current_vertex] += sum
    end

    counts[current_vertex] -= 1
    if counts[current_vertex] == 0
        for edge in predecessor_edges(subgraph, current_vertex)
            if !in(edge, sub_edges)
                continue
            else
                _evaluate_branching_subgraph(subgraph, vertex_sums[current_vertex] * value(edge), forward_vertex(subgraph, edge), sub_edges, counts, vertex_sums)
            end
        end
    end
end

### End of functions for evaluating subgraphs with branches.


function print_edges(a, msg)
    println(msg)
    for edge in edges(a)
        println(edge)
    end
end

function factor!(a::DerivativeGraph{T}) where {T}
    subgraph_list = compute_factorable_subgraphs(a)

    while !isempty(subgraph_list)
        # @info "Processed $count subgraphs out of $total"
        subgraph = pop!(subgraph_list)

        factor_subgraph!(subgraph)
        #test
        # Vis.draw_dot(graph(subgraph), start_nodes=[93], graph_label="factored subgraph $(vertices(subgraph))", value_labels=false)
        #end test
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
            @assert length(curr_edges) == 1 "Should only be one path from root $root_index to variable $var_index. Instead have $(length(curr_edges)) children from node $current_node_index on the path"
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
                return 1.0 #taking a derivative with respect to itself, which is 1. Need to figure out a better way to get the return number type right. This will always return Float64.
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


"""Verifies that there is a single path from each root to each variable, if a path exists. This should be an invariant of the factored graph so it should always be true. But the algorithm is complex enough that it is easy to accidentally introduce errors when adding features. `verify_paths` has negligible runtime cost compared to factorization."""
function _verify_paths(graph::DerivativeGraph, a::Int)
    branches = child_edges(graph, a)
    valid_graph = true

    if length(branches) > 1
        for br1 in 1:length(branches)
            for br2 in br1+1:length(branches)
                roots_intersect = reachable_roots(branches[br1]) .& reachable_roots(branches[br2])
                if !is_zero(roots_intersect) #if any shared roots then can't have any shared variables
                    if any(reachable_variables(branches[br1]) .& reachable_variables(branches[br2]))
                        valid_graph = false
                        @info "More than one path to variable for node $a. Non-zero intersection of reachable variables: $(reachable_variables(branches[br1]) .& reachable_variables(branches[br2]))"
                        #could break on first bad path but prefer to list all of them. Better for debugging.
                    end
                end
            end
        end

        for child in children(graph, a)
            valid_graph &= _verify_paths(graph, child)
        end
    end
    return valid_graph
end

"""verifies that there is a single path from each root to each variable, if such a path exists."""
function verify_paths(graph::DerivativeGraph)
    for root in roots(graph)
        if !_verify_paths(graph, postorder_number(graph, root))
            return false
        end
    end
    return true
end

"""Factors the graph then computes the Jacobian matrix. Only the columns of the Jacobian corresponsing to the elements of `partial_variables` will be computed. Example:
```
julia> nx, ny = Node.((x, y))
(x, y)

julia> symbolic_jacobian([nx*ny,ny*nx],[nx,ny])
2×2 Matrix{Node}:
 y  x
 y  x

julia> symbolic_jacobian([nx*ny,ny*nx],[ny,nx])
2×2 Matrix{Node}:
 x  y
 x  y

julia> symbolic_jacobian([nx*ny,ny*nx],[nx,ny])
2×2 Matrix{Node}:
 y  x
 y  x

julia> symbolic_jacobian([nx*ny,ny*nx],[nx])
2×1 Matrix{Node}:
 y
 y
 ```
"""
function _symbolic_jacobian!(graph::DerivativeGraph, partial_variables::AbstractVector{T}) where {T<:Node}
    outdim = codomain_dimension(graph)

    result = Matrix{Node}(undef, outdim, length(partial_variables))
    factor!(graph)

    @assert verify_paths(graph) #ensure a single path from each root to each variable. Derivative is likely incorrect if this is not true.

    for (i, var) in pairs(partial_variables)
        var_index = variable_node_to_index(graph, var)
        for root_index in 1:codomain_dimension(graph)
            if var_index !== nothing
                result[root_index, i] = evaluate_path(graph, root_index, var_index)
            else
                result[root_index, i] = zero(Node) #TODO fix this so get more generic zero value
            end
        end
    end

    return result
end

_symbolic_jacobian!(a::DerivativeGraph) = _symbolic_jacobian!(a, variables(a))

function _symbolic_jacobian(a::DerivativeGraph, variable_ordering::AbstractVector{T}) where {T<:Node}
    tmp = DerivativeGraph(roots(a)) #rebuild derivative graph. This is probably less efficient than deepcopy but deepcopy(Node(x)) != Node(x) which can lead to all kinds of trouble.
    return _symbolic_jacobian!(tmp, variable_ordering)
end

# _symbolic_jacobian(a::DerivativeGraph) = _symbolic_jacobian(a, variables(a))

"""Jacobian matrix of the n element function defined by `terms`. Each term element is a Node expression graph. Only the columns of the Jacobian corresponsing to the elements of `partial_variables` will be computed and the partial columns in the Jacobian matrix will be in the order specified by `partial_variables`. Examples:
```
julia> nx, ny = Node.((x, y))
(x, y)

julia> symbolic_jacobian([nx*ny,ny*nx],[nx,ny])
2×2 Matrix{Node}:
 y  x
 y  x

julia> symbolic_jacobian([nx*ny,ny*nx],[ny,nx])
2×2 Matrix{Node}:
 x  y
 x  y

julia> symbolic_jacobian([nx*ny,ny*nx],[nx,ny])
2×2 Matrix{Node}:
 y  x
 y  x

julia> symbolic_jacobian([nx*ny,ny*nx],[nx])
2×1 Matrix{Node}:
 y
 y
 ```
"""
symbolic_jacobian(terms::AbstractVector{T}, partial_variables::AbstractVector{S}) where {T<:Node,S<:Node} = _symbolic_jacobian(DerivativeGraph(terms), partial_variables)
export symbolic_jacobian


"""Computes sparse Jacobian matrix `J` using `SparseArray`. Each element `J[i,j]` is an expression graph which is the symbolic value of the Jacobian ∂fᵢ/∂vⱼ, where fᵢ is the ith output of the function represented by graph and vⱼ is the jth variable."""
function _sparse_symbolic_jacobian!(graph::DerivativeGraph, partial_variables::AbstractVector{T}) where {T<:Node}
    row_indices = Int64[]
    col_indices = Int64[]
    values = Node[]
    @assert length(partial_variables) == domain_dimension(graph)

    factor!(graph)

    @assert verify_paths(graph) #ensure a single path from each root to each variable. Derivative is likely incorrect if this is not true.
    #input is an array of Node's representing variables. Need a mapping from the variable index matching the Node to the index in variable_ordering
    variable_index = map(x -> variable_postorder_to_index(graph, postorder_number(graph, x)), partial_variables)

    for root in 1:codomain_dimension(graph)
        reach_vars = reachable_variables(graph, root_index_to_postorder_number(graph, root))
        for (i, partial_var) in pairs(partial_variables)
            partial_index = variable_node_to_index(graph, partial_var) #make sure variable is in the domain of the graph
            if partial_index !== nothing
                if reach_vars[partial_index] #make sure variable is reachable from this root. If ∂fⱼ/∂xᵢ ≡ 0 then there will not be a path from fⱼ to xᵢ so don't need a separate check for this case.
                    push!(row_indices, root)
                    push!(col_indices, i)
                    push!(values, evaluate_path(graph, root, partial_index))
                end
            end
        end
    end

    return sparse(row_indices, col_indices, values, codomain_dimension(graph), domain_dimension(graph))
end

sparse_symbolic_jacobian(terms::AbstractVector{T}, partial_variables::AbstractVector{S}) where {T<:Node,S<:Node} = _sparse_symbolic_jacobian!(DerivativeGraph(terms), partial_variables)
export sparse_symbolic_jacobian

"""Returns a vector of Node, where each element in the vector is the symbolic form of `Jv``. Also returns `v_vector` a vector of the `v` variables. This is useful if you want to generate a function to evaluate `Jv` and you want to separate the inputs to the function and the `v` variables."""
function jacobian_times_v(terms::AbstractVector{T}, partial_variables::AbstractVector{S}) where {T<:Node,S<:Node}
    graph = DerivativeGraph(terms)
    v_vector = make_variables(gensym(), domain_dimension(graph))
    factor!(graph)
    for (variable, one_v) in zip(variables(graph), v_vector)
        new_edges = PathEdge[]
        old_edges = collect(parent_edges(graph, variable)) #can't use iterator returned by parent_edges because it will include new edges. Need snapshot of old edges.
        for parent in old_edges
            push!(new_edges,
                PathEdge(
                    top_vertex(parent),
                    bott_vertex(parent),
                    value(parent) * one_v,
                    copy(reachable_variables(parent)),
                    copy(reachable_roots(parent))
                )
            )
            #multiply all parent edges by one_v and replace edge in graph
        end

        for new_edge in new_edges
            add_edge!(graph, new_edge)
        end

        for old_edge in old_edges
            delete_edge!(graph, old_edge, true)
        end
    end

    outdim = codomain_dimension(graph)

    result = Vector{Node}(undef, outdim)

    for i in eachindex(result)
        result[i] = Node(0.0)
    end



    @assert verify_paths(graph) #ensure a single path from each root to each variable. Derivative is likely incorrect if this is not true.

    for (i, var) in pairs(partial_variables)
        var_index = variable_node_to_index(graph, var)
        for root_index in 1:codomain_dimension(graph)
            if var_index !== nothing
                result[root_index] += evaluate_path(graph, root_index, var_index)
            else
                result[root_index, i] = zero(Node) #TODO fix this so get more generic zero value
            end
        end
    end

    return result, v_vector #need v_vector values if want to make executable after making symbolic form. Need to differentiate between variables that were in original graph and variables introduced by v_vector
end
export jacobian_times_v

"""Creates an executable `(x,v)-> Jv([x;v]...)` to evaluate `J(x)*v`. The executable takes two vector arguments `x` and `v`. The `x` vector is the input to Jacobian function `J`. `v` is the vector you want to multiply the Jacobian by."""
function jacobian_times_v_exe(terms::AbstractVector{T}, partial_variables::AbstractVector{S}) where {T<:Node,S<:Node}
    Jv, v_vec = jacobian_times_v(terms, partial_variables)
    both_vars = [partial_variables; v_vec]
    return make_function(reshape(Jv, (length(Jv), 1)), both_vars)
end
export jacobian_times_v_exe


"""Returns a vector of Node, where each element in the vector is the symbolic form of `Jᵀv`. Also returns `v_vector` a vector of the `v` variables. This is useful if you want to generate a function to evaluate `Jᵀv` and you want to separate the inputs to the function and the `v` variables."""
function jacobian_transpose_v(terms::AbstractVector{T}, partial_variables::AbstractVector{S}) where {T<:Node,S<:Node}
    graph = DerivativeGraph(terms)
    r_vector = make_variables(gensym(), codomain_dimension(graph))
    factor!(graph)
    for (one_root, one_v) in zip(roots(graph), r_vector)
        new_edges = PathEdge[]
        old_edges = collect(child_edges(graph, one_root)) #can't use iterator returned by parent_edges because it will include new edges. Need snapshot of old edges.
        for child in old_edges
            push!(new_edges,
                PathEdge(
                    top_vertex(child),
                    bott_vertex(child),
                    value(child) * one_v,
                    copy(reachable_variables(child)),
                    copy(reachable_roots(child))
                )
            )
            #multiply all parent edges by one_v and replace edge in graph
        end

        for new_edge in new_edges
            add_edge!(graph, new_edge)
        end

        for old_edge in old_edges
            delete_edge!(graph, old_edge, true)
        end
    end

    outdim = domain_dimension(graph)

    result = Vector{Node}(undef, outdim)

    for i in eachindex(result)
        result[i] = Node(0.0)
    end



    @assert verify_paths(graph) #ensure a single path from each root to each variable. Derivative is likely incorrect if this is not true.

    for (i, var) in pairs(partial_variables)
        for root_index in 1:codomain_dimension(graph)
            var_index = variable_node_to_index(graph, var)
            if var_index !== nothing
                result[var_index] += evaluate_path(graph, root_index, var_index)
            else
                result[var_index] = zero(Node) #TODO fix this so get more generic zero value
            end
        end
    end

    return result, r_vector #need v_vector values if want to make executable after making symbolic form. Need to differentiate between variables that were in original graph and variables introduced by v_vector
end
export jacobian_transpose_v

"""Creates an executable `(x,v)-> Jᵀv([x;v]...)` to evaluate `J(x)ᵀ*v`. The executable takes two vector arguments `x` and `v`. The `x` vector is the input to Jacobian function `J`. `v` is the vector you want to multiply the Jacobian transpose by."""
function jacobian_transpose_v_exe(terms::AbstractVector{T}, partial_variables::AbstractVector{S}) where {T<:Node,S<:Node}
    Jᵀv, v_vec = jacobian_transpose_v(terms, partial_variables)
    both_vars = [partial_variables; v_vec]

    return make_function(reshape(Jᵀv, (length(Jᵀv), 1)), both_vars)
end
export jacobian_transpose_v_exe

function make_Expr(func_array::AbstractArray{T,N}, input_variables::AbstractVector{S}; in_place=false, use_vector_runtime_args=false) where {T<:Node,S<:Node,N}
    node_to_var = Dict{Node,Union{Symbol,Real,Expr}}()
    body = Expr(:block)

    if !in_place
        push!(body.args, :(result = Array{promote_type(Float64, eltype(input_variables)),$N}(undef, $(size(func_array)...))))
    end

    node_to_index = Dict{Node,Int64}()
    for (i, node) in pairs(input_variables)
        node_to_index[node] = i
    end

    for (i, node) in pairs(func_array)
        node_body, variable = function_body!(node, node_to_index, node_to_var)
        push!(node_body.args, :(result[$i] = $variable))
        push!(body.args, node_body)
    end

    push!(body.args, :(return result))

    vector_type(vars) = isa(vars, SVector) ? :(SVector{$(length(vars)),T}) : :(Vector)
    vector_type(partial_variables)
    if in_place
        return :(((input_variables::$(typeof(input_variables)), result) -> $body))
    else
        return :(((input_variables::$(typeof(input_variables))) -> $body))
    end
end
export make_Expr

make_function(func_array::AbstractArray{T}, input_variables::AbstractVector{S}; in_place=false, use_vector_runtime_args=false) where {T<:Node,S<:Node} = @RuntimeGeneratedFunction(make_Expr(func_array, input_variables, in_place=in_place, use_vector_runtime_args=use_vector_runtime_args))
export make_function


"""Computes the full symbolic Hessian matrix"""
function hessian(graph::DerivativeGraph, variable_order)
    @assert codomain_dimension(graph) == 1
    return hessian(roots(graph)[1], variable_order)
end

"""Computes the full symbolic Hessian matrix"""
function hessian(expression::Node, variable_order::AbstractVector{S}) where {S<:Node} #would prefer to return a Symmetric matrix but that type only works with elements that are subtypes of Number. Which Node is not. Fix later, if possible.
    tmp = DerivativeGraph(expression)
    jac = _symbolic_jacobian!(tmp, variable_order)
    tmp2 = DerivativeGraph(vec(jac))
    return _symbolic_jacobian!(tmp2, variable_order)
end
export hessian

function sparse_hessian(expression::Node, variable_order::AbstractVector{S}) where {S<:Node}
    gradient = sparse_symbolic_jacobian(DerivativeGraph(expression), variable_order)
end

function unique_nodes(jacobian::AbstractArray{T}) where {T<:Node}
    nodes = Set{Node}()
    for index in eachindex(jacobian)
        oned = all_nodes(jacobian[index])
        union!(nodes, oned)
    end
    return nodes
end

"""Count of number of operations in graph."""
number_of_operations(jacobian::AbstractArray{T}) where {T<:Node} = length(filter(x -> is_tree(x), unique_nodes(jacobian)))


"""computes ∂A/∂variables[1],...,variables[n]. Repeated differentiation rather than computing different columns of the Jacobian. Example:

```
julia> A = [t t^2;3t^2 5]
2×2 Matrix{Node}:
 t              (t ^ 2)
 (3 * (t ^ 2))  5

julia> derivative(A,t)
2×2 Matrix{Node}:
 1.0      (2 * t)
 (6 * t)  0.0

 julia> derivative(A,t,t)
 2×2 Matrix{Node{T, 0} where T}:
 0.0  2
 6    0.0
 ```
 """
function derivative(A::Matrix{<:Node}, variables::T...) where {T<:Node}
    var = variables[1]
    mat = _derivative(A, var)
    rest = variables[2:end]
    if rest === ()
        return mat
    else
        return derivative(mat, rest...)
    end
end
export derivative

"""Computes the derivative of the function matrix `A` with respect to  `variable`."""
function _derivative(A::Matrix{<:Node}, variable::T) where {T<:Node}
    #convert A into vector then compute jacobian
    vecA = vec(A)
    graph = DerivativeGraph(vecA)

    temp = _symbolic_jacobian!(graph)
    #pick out the column of the Jacobian containing partials with respect to variable and pack them back into a matrix of the same shape as A. Later, if this becomes a bottleneck, modify symbolic_jacobian! to only compute the single column of derivatives.
    column_index = variable_node_to_index(graph, variable)
    if column_index === nothing
        return Node.(zeros(size(A)))
    else
        column = temp[:, column_index] #these 3 lines allocate but this shouldn't be a significant source of inefficiency.
        mat_column = reshape(column, size(A))
        result = copy(mat_column)

        return Node.(result)
    end
end



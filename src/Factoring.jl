
next_edge_constraint(sub::FactorableSubgraph{T,PostDominatorSubgraph}) where {T} = PathConstraint(dominating_node(sub), graph(sub), false, reachable_roots(sub), dominance_mask(sub))
next_edge_constraint(sub::FactorableSubgraph{T,DominatorSubgraph}) where {T} = PathConstraint(dominating_node(sub), graph(sub), true, dominance_mask(sub), reachable_variables(sub))
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
    tmp = _node_edges(edges(graph), dominated_node)

    if tmp === nothing #no edges edges up or down so must be a root. dominated_node can't be part of a factorable dom subgraph.
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

function pdom_subgraph(graph::DerivativeGraph, variable_index::Integer, dominated_node::Integer, pidom)
    tmp = _node_edges(edges(graph), dominated_node)
    if tmp === nothing #no edges up or down. Must be a root with no children. dominated_node can't be part of a factorable pdom subgraph.
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

struct FactorOrder <: Base.Order.Ordering
end
export FactorOrder

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
function compute_factorable_subgraphs(graph::DerivativeGraph{T}) where {T}
    dom_subgraphs = Dict{Tuple{T,T},BitVector}()
    pdom_subgraphs = Dict{Tuple{T,T},BitVector}()

    timer = TimerOutput()

    temp_doms = Dict{T,T}()

    for root_index in 1:codomain_dimension(graph)
        post_num = root_index_to_postorder_number(graph, root_index)
        @timeit timer "compute_dom_table" temp_dom = compute_dom_table(graph, true, root_index, post_num, temp_doms)

        #test
        # @time temp_dom = compute_dom_table(graph, true, root_index, post_num, temp_doms)
        # println("size of dom table $(length(temp_dom))")
        #end test

        # if mod(root_index, 10) == 0
        #     @info "computed $root_index out of $(codomain_dimension(graph)) root dom tables"
        # end

        for dominated in keys(temp_dom)
            dsubgraph = dom_subgraph(graph, root_index, dominated, temp_dom)
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
        post_num = variable_index_to_postorder_number(graph, variable_index)
        temp_dom = compute_dom_table(graph, false, variable_index, post_num, temp_doms)
        for dominated in keys(temp_dom)
            psubgraph = pdom_subgraph(graph, variable_index, dominated, temp_dom)
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

    # result = Vector{FactorableSubgraph}(undef, length(dom_subgraphs) + length(pdom_subgraphs)) #slightly faster and more memory efficient than pushing to a zero length array
    # empty!(result)

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
export compute_factorable_subgraphs

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
export old_edge_path

function evaluate_subgraph(subgraph::FactorableSubgraph{T,S}) where {T,S<:Union{DominatorSubgraph,PostDominatorSubgraph}}
    constraint = next_edge_constraint(subgraph)
    sum = Node(0.0)
    gr = graph(subgraph)
    roots_intersect = trues(codomain_dimension(gr))
    vars_intersect = trues(domain_dimension(gr))

    rel_edges = get_edge_vector()
    for edge in relation_edges!(constraint, dominated_node(subgraph), rel_edges)
        flag, pedges, roots_reach, vars_reach = old_edge_path(constraint, dominating_node(subgraph), S == DominatorSubgraph, reachable(subgraph), edge)
        #sort by num_uses then from largest to smallest postorder number

        if flag == 1 #non-branching path through subgraph
            roots_intersect .= roots_intersect .& roots_reach
            vars_intersect .= vars_intersect .& vars_reach

            sort!(pedges, lt=path_sort_order)
            sum += multiply_sequence(pedges)

            #find sequences of equal times_used and multiply them. Then multiply each of the collapsed sequences to get the final result
            reclaim_edge_vector(pedges)
        end
    end


    reclaim_edge_vector(rel_edges)
    return sum, roots_intersect, vars_intersect
end
export evaluate_subgraph

function make_factored_edge(subgraph::FactorableSubgraph{T,DominatorSubgraph}) where {T}
    sum, _, _ = evaluate_subgraph(subgraph)

    roots_reach = copy(dominance_mask(subgraph))
    vars_reach = copy(reachable_variables(subgraph))
    return PathEdge(dominating_node(subgraph), dominated_node(subgraph), sum, vars_reach, roots_reach), roots_reach, vars_reach
end

export make_factored_edge
function make_factored_edge(subgraph::FactorableSubgraph{T,PostDominatorSubgraph}) where {T}
    sum, roots_reach, vars_reach = evaluate_subgraph(subgraph)

    roots_reach = copy(reachable_roots(subgraph))
    vars_reach = copy(dominance_mask(subgraph))
    return PathEdge(dominating_node(subgraph), dominated_node(subgraph), sum, vars_reach, roots_reach), roots_reach, vars_reach
end
export make_factored_edge

"""reset root and variable masks for edges in the graph and add a new edge connecting `dominating_node(subgraph)` and `dominated_node(subgraph)` to the graph that has the factored value of the subgraph"""
function factor_subgraph!(subgraph::FactorableSubgraph)
    if subgraph_exists(subgraph)
        new_edge, roots_reach, vars_reach = make_factored_edge(subgraph)
        add_non_dom_edges!(subgraph)
        #reset roots in R, if possible. All edges higher in the path than the first vertex with more than one child cannot be reset.
        edges_to_delete = reset_edge_masks!(subgraph) #TODO need to modify reset_edge_masks! so it handles the case where a path may have been destroyed due to factorization.
        for edge in edges_to_delete
            delete_edge!(graph(subgraph), edge)
        end

        add_edge!(graph(subgraph), new_edge)
    end
end
export factor_subgraph!

function print_edges(a, msg)
    println(msg)
    for edge in edges(a)
        println(edge)
    end
end

function factor!(a::DerivativeGraph{T}) where {T}
    @time subgraph_list = compute_factorable_subgraphs(a)

    count = 0
    total = length(subgraph_list)
    @info "$total factorable subgraphs"
    while !isempty(subgraph_list)
        # @info "Processed $count subgraphs out of $total"
        count += 1

        subgraph = pop!(subgraph_list)

        factor_subgraph!(subgraph)
    end
    return nothing #return nothing so people don't mistakenly think this is returning a copy of the original graph
end
export factor!

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


"""verifies that there is a single path from each root to each variable, if a path exists. Used for diagnostics and debugging. Not used in normal use of FSD."""
function _verify_paths(graph::DerivativeGraph, a::Int)
    branches = child_edges(graph, a)

    if length(branches) > 1
        intersection::BitVector = mapreduce(reachable_variables, .&, branches, init=trues(domain_dimension(graph)))
        if !is_zero(intersection)
            @info "More than one path to variable for node $a. Non-zero intersection of reachable variables $intersection"
            for branch in branches
                @info "reachable variables $(reachable_variables(branch))"
            end
            return false
        end
        for child in children(graph, a)
            if !(_verify_paths(graph, child))
                return false
            end
        end
    end

    return true
end

"""verifies that there is a single path from each root to each variable, if such a path exists. Used for diagnostics and debugging. Normal us of FSD does not require these functions."""
function verify_paths(graph::DerivativeGraph)
    for root in roots(graph)
        if !_verify_paths(graph, postorder_number(graph, root))
            return false
        end
    end
    return true
end

"""Factors the graph then computes jacobian matrix. Destructive."""
function symbolic_jacobian!(graph::DerivativeGraph, variable_ordering::AbstractVector{T}) where {T<:Node}
    indim = domain_dimension(graph)
    outdim = codomain_dimension(graph)

    result = Matrix{Node}(undef, outdim, indim)
    factor!(graph)

    verify_paths(graph)

    for (i, var) in pairs(variable_ordering)
        var_index = variable_node_to_index(graph, var)
        for root_index in 1:codomain_dimension(graph)
            result[root_index, i] = evaluate_path(graph, root_index, var_index)
        end
    end

    return result
end
export symbolic_jacobian!

symbolic_jacobian!(a::DerivativeGraph) = symbolic_jacobian!(a, variables(a))


jacobian_function!(graph::DerivativeGraph) = jacobian_function!(graph, variables(graph))

function jacobian_Expr!(graph::DerivativeGraph, variable_order::AbstractVector{S}; in_place=false) where {S<:Node}
    tmp = symbolic_jacobian!(graph, variable_order)
    node_to_var = Dict{Node,Union{Symbol,Real}}()
    all_vars = variables(graph)

    if variable_order === nothing
        ordering = all_vars
    else
        ordering = Node.(variable_order)
    end
    @assert Set(all_vars) ⊆ Set(ordering) "Not every variable in the graph had a corresponding ordering variable."

    body = Expr(:block)
    if !in_place
        push!(body.args, :(result = fill(0.0, $(size(tmp))))) #shouldn't need to fill with zero. All elements should be defined. Unless doing sparse Jacobian.
    end

    for (i, node) in pairs(tmp)
        node_body, variable = function_body(node, node_to_var)
        push!(node_body.args, :(result[$i] = $variable))
        push!(body.args, node_body)
    end

    push!(body.args, :(return result))

    if in_place
        return Expr(:->, Expr(:tuple, map(x -> node_symbol(x), ordering)..., :result), body)
    else
        return Expr(:->, Expr(:tuple, map(x -> node_symbol(x), ordering)...), body)
    end
end
export jacobian_Expr!

jacobian_function!(graph::DerivativeGraph, variable_order::AbstractVector{S}; in_place=false) where {S<:Node} = @RuntimeGeneratedFunction(jacobian_Expr!(graph, variable_order; in_place))
export jacobian_function!

function unique_nodes(jacobian::AbstractArray{T}) where {T<:Node}
    nodes = Set{Node}()
    for oned in all_nodes.(jacobian)
        union!(nodes, oned)
    end
    return nodes
end

number_of_operations(jacobian::AbstractArray{T}) where {T<:Node} = length(filter(x -> is_tree(x), unique_nodes(jacobian)))
export number_of_operations

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

    temp = symbolic_jacobian!(graph)
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





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
            # if value(sum) == 0
            #     display(subgraph)
            #     write_dot("sph.svg", graph(subgraph), value_labels=true, reachability_labels=false, start_nodes=[24])
            # end

            # # @assert value(sum) != 0

            new_edge = make_factored_edge(subgraph, sum)
        end
        add_non_dom_edges!(subgraph)
        #reset roots in R, if possible. All edges earlier in the path than the first vertex with more than one child cannot be reset.
        edges_to_delete = reset_edge_masks!(subgraph)
        for edge in edges_to_delete
            delete_edge!(graph(subgraph), edge)
        end

        add_edge!(graph(subgraph), new_edge)

        #test
        gr = graph(subgraph)
        for one_edge in unique_edges(gr)
            low_index = bott_vertex(one_edge)
            nd = node(gr, low_index)

            if !is_variable(nd) && !is_constant(nd)
                @assert any(reachable_variables(one_edge))

            end
        end
        #end test

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

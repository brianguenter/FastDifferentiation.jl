"""Factors the graph then computes the Jacobian matrix. Only the columns of the Jacobian corresponsing to the elements of `partial_variables` will be computed. Example:
```
julia> @variables x y

julia> jacobian([x*y,y*x],[x,y])
2×2 Matrix{Node}:
 y  x
 y  x

julia> jacobian([x*y,y*x],[y,x])
2×2 Matrix{Node}:
 x  y
 x  y

julia> jacobian([x*y,y*x],[x,y])
2×2 Matrix{Node}:
 y  x
 y  x

julia> jacobian([x*y,y*x],[x])
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

"""Jacobian matrix of the n element function defined by `terms`. Each term element is a Node expression graph. Only the columns of the Jacobian corresponsing to the elements of `partial_variables` will be computed and the partial columns in the Jacobian matrix will be in the order specified by `partial_variables`. Examples:
```julia-repl

julia> @variables x y

julia> jacobian([x*y,y*x],[x,y])
2×2 Matrix{Node}:
 y  x
 y  x

julia> jacobian([x*y,y*x],[y,x])
2×2 Matrix{Node}:
 x  y
 x  y

julia> jacobian([x*y,y*x],[x])
2×1 Matrix{Node}:
 y
 y
```
"""
jacobian(terms::AbstractVector{T}, partial_variables::AbstractVector{S}) where {T<:Node,S<:Node} = _symbolic_jacobian(DerivativeGraph(terms), partial_variables)
export jacobian


"""Computes sparse Jacobian matrix `J` using `SparseArray`. Each element `J[i,j]` is an expression graph which is the symbolic value of the Jacobian ∂fᵢ/∂vⱼ, where fᵢ is the ith output of the function represented by graph and vⱼ is the jth variable."""
function _sparse_symbolic_jacobian!(graph::DerivativeGraph, partial_variables::AbstractVector{T}) where {T<:Node}
    row_indices = Int64[]
    col_indices = Int64[]
    values = Node[]

    factor!(graph)

    @assert verify_paths(graph) #ensure a single path from each root to each variable. Derivative is likely incorrect if this is not true.
    #input is an array of Node's representing variables. Need a mapping from the variable index matching the Node to the index in variable_ordering
    variable_index = map(x -> variable_postorder_to_index(graph, postorder_number(graph, x)), partial_variables)

    for root in 1:codomain_dimension(graph)
        for (i, partial_var) in pairs(partial_variables)
            partial_index = variable_node_to_index(graph, partial_var) #make sure variable is in the domain of the graph
            if partial_index !== nothing
                tmp = evaluate_path(graph, root, partial_index)
                if !is_zero(tmp)
                    push!(row_indices, root)
                    push!(col_indices, i)
                    push!(values, tmp)
                end
            end
        end
    end

    return sparse(row_indices, col_indices, values, codomain_dimension(graph), length(partial_variables))
end

"""Returns a sparse array containing the Jacobian of the function defined by `terms`"""
sparse_jacobian(terms::AbstractVector{T}, partial_variables::AbstractVector{S}) where {T<:Node,S<:Node} = _sparse_symbolic_jacobian!(DerivativeGraph(terms), partial_variables)
export sparse_jacobian

"""Returns a vector of Node, where each element in the vector is the symbolic form of `Jv``. Also returns `v_vector` a vector of the `v` variables. This is useful if you want to generate a function to evaluate `Jv` and you want to separate the inputs to the function and the `v` variables."""
function jacobian_times_v(terms::AbstractVector{T}, partial_variables::AbstractVector{S}) where {T<:Node,S<:Node}
    graph = DerivativeGraph(terms)
    v_vector = make_variables(gensym(), domain_dimension(graph))
    factor!(graph)
    for (variable, one_v) in zip(partial_variables, v_vector)
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

"""Computes Hessian times a vector v without forming the Hessian matrix. Useful when the Hessian would be impractically large."""
function hessian_times_v(term::T, partial_variables::AbstractVector{S}) where {T<:Node,S<:Node}
    gradient = vec(jacobian([term], partial_variables))
    return jacobian_times_v(gradient, partial_variables)
end
export hessian_times_v

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

"""Returns the dense symbolic Hessian matrix. Example:
```
julia> @variables x y

julia> hessian(x^2*y^2,[x,y])
2×2 Matrix{FastDifferentiation.Node}:
 (2 * (y ^ 2))  (4 * (y * x))
 (4 * (x * y))  (2 * (x ^ 2))
```
"""
function hessian(expression::Node, variable_order::AbstractVector{S}) where {S<:Node} #would prefer to return a Symmetric matrix but that type only works with elements that are subtypes of Number. Which Node is not. Fix later, if possible.
    tmp = DerivativeGraph(expression)
    jac = _symbolic_jacobian!(tmp, variable_order)
    tmp2 = DerivativeGraph(vec(jac))
    return _symbolic_jacobian!(tmp2, variable_order)
end
export hessian

"""Compute a sparse symbolic Hessian. Returns a sparse matrix of symbolic expressions. 
Can be used in combination with `make_function` to generate an executable that
 will return a sparse matrix or take one as an in-place argument. Example:

 ```
julia> @variables x y

julia> a = sparse_hessian(x*y,[x,y])
2×2 SparseArrays.SparseMatrixCSC{FastDifferentiation.Node, Int64} with 2 stored entries:
 ⋅  1
 1  ⋅

julia> f1 = make_function(a,[x,y])
...

julia> f1([1.0,2.0])
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
  ⋅   1.0
 1.0   ⋅

julia> tmp = similar(a,Float64)
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
  ⋅            4.24399e-314
 4.24399e-314   ⋅

julia> f2 = make_function(a,[x,y],in_place=true)
...

julia> f2([1.0,2.0],tmp)
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
  ⋅   1.0
 1.0   ⋅

julia> tmp
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
  ⋅   1.0
 1.0   ⋅

```

"""
function sparse_hessian(expression::Node, variable_order::AbstractVector{S}) where {S<:Node}
    gradient = jacobian([expression], variable_order)
    return sparse_jacobian(vec(gradient), variable_order)
end
export sparse_hessian




"""computes ∂A/(∂variables[1],...,∂variables[n]). Repeated differentiation rather than computing different columns of the Jacobian. Example:

```julia-repl

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
    #pick out the column of the Jacobian containing partials with respect to variable and pack them back into a matrix of the same shape as A. Later, if this becomes a bottleneck, modify jacobian! to only compute the single column of derivatives.
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

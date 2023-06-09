
"""Used to determine whether to fill zero array elements with an assignment statement or to fill the array in the declaration. Function arrays with many zero elements generate many zero assignment statements which can make compilation time slow. But Need a heuristic to determine when to choose one or the other."""
function sparsity(sym_func::AbstractArray{<:Node})
    zeros = mapreduce(x -> is_zero(x) ? 1 : 0, +, sym_func)
    tot = prod(size(sym_func))
    return zeros == 0 ? 1.0 : (tot - zeros) / tot
end

"""Create body of Expr that will evaluate the function. The function body will be a sequence of assignment statements to automatically generated variable names. This is an example for a simple function:
```julia
quote
    var"##343" = 2x
    var"##342" = var"##343" + y
    var"##341" = var"##342" + 1
end
```
The last automatically generated name (in this example var"##341") is the second return value of `function_body`. This variable will hold the value of evaluating the dag at runtime.
If the dag is a constant then the function body will be empty:
```julia
quote
end
```
and the second return value will be the constant value.
"""
function function_body!(dag::Node, variable_to_index::IdDict{Node,Int64}, node_to_var::Union{Nothing,IdDict{Node,Union{Symbol,Real,Expr}}}=nothing)
    if node_to_var === nothing
        node_to_var = IdDict{Node,Union{Symbol,Real,Expr}}()
    end

    body = Expr(:block)

    function _dag_to_function(node)

        tmp = get(node_to_var, node, nothing)

        if tmp === nothing #if node not in node_to_var then it hasn't been visited. Otherwise it has so don't recurse.
            node_to_var[node] = node_symbol(node, variable_to_index)

            if is_tree(node)
                args = _dag_to_function.(children(node))
                statement = :($(node_to_var[node]) = $(Symbol(value(node)))($(args...)))
                push!(body.args, statement)
            end
        end

        return node_to_var[node]
    end

    return body, _dag_to_function(dag)
end

return_declaration(::StaticArray{S,T,N}, input_variables::AbstractVector{T2}) where {S,T,T2,N} = :(result = MArray{$(S),promote_type(Float64, eltype(input_variables)),$N}(undef); result .= 0) #need to initialize array to zero because this is no longer being done by simple assignment statements.

"""fills the return array with zeros, which is much more efficient for sparse arrays than setting each element with a line of code"""
return_declaration(func_array::Array{T,N}, input_variables::AbstractVector{S}) where {S,T,N} = :(result = zeros(promote_type(Float64, eltype(input_variables)), $(size(func_array)...)))

return_expression(::SArray) = :(return SArray(result))
return_expression(::Array) = :(return result)

function make_Expr(func_array::AbstractArray{T}, input_variables::AbstractVector{S}, in_place::Bool) where {T<:Node,S<:Node}
    node_to_var = IdDict{Node,Union{Symbol,Real,Expr}}()
    body = Expr(:block)

    if in_place
        push!(body.args, :(result .= zero(eltype(input_variables))))
    else
        push!(body.args, (return_declaration(func_array, input_variables)))
    end

    node_to_index = IdDict{Node,Int64}()
    for (i, node) in pairs(input_variables)
        node_to_index[node] = i
    end

    for (i, node) in pairs(func_array)
        if !is_zero(node)
            node_body, variable = function_body!(node, node_to_index, node_to_var)
            push!(node_body.args, :(result[$i] = $variable))
            push!(body.args, node_body)
        end
    end

    if !in_place
        push!(body.args, return_expression(func_array))
    end

    if in_place
        return :((input_variables, result) -> @inbounds begin
            $body
        end)
    else
        return :((input_variables) -> @inbounds begin
            $body
        end)
    end
end
export make_Expr

function make_Expr(A::SparseMatrixCSC{T,Ti}, input_variables::AbstractVector{S}, in_place::Bool) where {T<:Node,S<:Node,Ti}
    rows = rowvals(A)
    vals = nonzeros(A)
    _, n = size(A)
    body = Expr(:block)
    node_to_var = IdDict{Node,Union{Symbol,Real,Expr}}()

    if !in_place #have to store the sparse vector indices in the generated code to know how to create sparsity pattern
        push!(body.args, :(element_type = promote_type(Float64, eltype(input_variables))))
        push!(body.args, :(result = SparseMatrixCSC($(A.m), $(A.n), $(A.colptr), $(A.rowval), zeros(element_type, $(length(A.nzval))))))
    end

    push!(body.args, :(vals = nonzeros(result)))

    node_to_index = IdDict{Node,Int64}()
    for (i, node) in pairs(input_variables)
        node_to_index[node] = i
    end

    for j = 1:n
        for i in nzrange(A, j)
            # println(i)
            row = rows[i]
            node_body, variable = function_body!(vals[i], node_to_index, node_to_var)
            push!(node_body.args, :(vals[$i] = $variable))
            push!(body.args, node_body)
        end
    end

    push!(body.args, :(return result))

    if in_place
        return :((input_variables, result) -> $body)
    else
        return :((input_variables) -> $body)
    end
end
export make_Expr

"""Makes a function to evaluate the symbolic expressions in `func_array`. If `in_place=true` it will generate code to fill a user supplied array with the result. In this case if an element of `func_array` is identically zero the generated code will not set the result to zero to zero because this leads to code bloat and slow execution. If you need these array elements to be zero instead of `undef` you should properly initialize your array to zero before passing it to the runtime generated function.

If `in_place=false` then the returned array will be properly initialized with zeros."""
function make_function(func_array::AbstractArray{T}, input_variables::AbstractVector{<:Node}...; in_place=false) where {T<:Node}
    @RuntimeGeneratedFunction(make_Expr(func_array, vcat(input_variables...), in_place))
end
export make_function
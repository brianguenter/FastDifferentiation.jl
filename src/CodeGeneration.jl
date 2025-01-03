
"""
    sparsity(sym_func::AbstractArray{<:Node})


Computes a number representing the sparsity of the array of expressions. If `nelts` is the number of elements in the array and `nzeros` is the number of zero elements in the array
then `sparsity = (nelts-nzeros)/nelts`. 

Frequently used in combination with a call to `make_function` to determine whether to set keyword argument `init_with_zeros` to false."""
function sparsity(sym_func::AbstractArray{<:Node})
    zeros = mapreduce(x -> is_identically_zero(x) ? 1 : 0, +, sym_func)
    tot = prod(size(sym_func))
    return zeros == 0 ? 1.0 : (tot - zeros) / tot
end
export sparsity


function _dag_to_function!(node, local_body, variable_to_index, node_to_var)

    tmp = get(node_to_var, node, nothing)

    if tmp === nothing #if node not in node_to_var then it hasn't been visited. Otherwise it has so don't recurse.
        node_to_var[node] = node_symbol(node, variable_to_index)

        if is_tree(node)
            if value(node) === if_else #special case code generation for if...else. Need to generate nested code so only the statements in the true or false branch will be executed.
                true_body = Expr(:block)
                false_body = Expr(:block)
                if_cond_var = _dag_to_function!(children(node)[1], local_body, variable_to_index, node_to_var)

                true_node = children(node)[2]
                false_node = children(node)[3]

                if is_leaf(true_node) #handle leaf nodes properly 
                    if is_constant(true_node)
                        temp_val = value(true_node)
                    else
                        temp_val = node_to_var[true_node]
                    end

                    push!(true_body.args, :($(gensym(:s)) = $(temp_val))) #seems roundabout to use an assignment when really just want the value of the node but couldn't figure out how to make this work with Expr
                else
                    _dag_to_function!(children(node)[2], true_body, variable_to_index, node_to_var)
                end

                if is_leaf(false_node)
                    if is_constant(false_node)
                        temp_val = value(false_node)
                    else
                        temp_val = node_to_var[false_node]
                    end
                    push!(false_body.args, :($(gensym(:s)) = $(temp_val))) #seems roundabout to use an assignment when really just want the value of the node but couldn't figure out how to make this work with Expr
                else
                    _dag_to_function!(children(node)[3], false_body, variable_to_index, node_to_var)
                end

                statement = :($(node_to_var[node]) = if $(if_cond_var)
                    $(true_body)
                else
                    $(false_body)
                end)
            else
                args = _dag_to_function!.(children(node), Ref(local_body), Ref(variable_to_index), Ref(node_to_var))
                statement = :($(node_to_var[node]) = $(Symbol(value(node)))($(args...)))
            end
            push!(local_body.args, statement)
        end
    end

    return node_to_var[node]
end

"""
    function_body!(
        dag::Node,
        variable_to_index::IdDict{Node,Int64},
        node_to_var::Union{Nothing,IdDict{Node,Union{Symbol,Real,Expr}}}=nothing
    )

Create body of Expr that will evaluate the function. The function body will be a sequence of assignment statements to automatically generated variable names. This is an example for a simple function:
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
function function_body!(dag::Node, variable_to_index::IdDict{Node,Union{Expr,Int64}}, node_to_var::Union{Nothing,IdDict{Node,Union{Symbol,Real,Expr}}}=nothing, body::Expr=Expr(:block))
    if node_to_var === nothing
        node_to_var = IdDict{Node,Union{Symbol,Real,Expr}}()
    end

    return body, _dag_to_function!(dag, body, variable_to_index, node_to_var)
end

function zero_array_declaration(array::StaticArray{S,<:Any,N}) where {S,N}
    #need to initialize array to zero because this is no longer being done by simple assignment statements.
    :($(undef_array_declaration(array)); result .= 0)
end

function undef_array_declaration(::StaticArray{S,<:Any,N}) where {S,N}
    #need to initialize array to zero because this is no longer being done by simple assignment statements.
    :(result = MArray{$(S),result_element_type,$N}(undef))
end

"""
    return_declaration(func_array::Array, input_variables::AbstractVector)

Fills the return array with zeros, which is much more efficient for sparse arrays than setting each element with a line of code
"""
zero_array_declaration(func_array::Array{T,N}) where {T,N} = :(result = zeros(result_element_type, $(size(func_array))))
undef_array_declaration(func_array::Array{T,N}) where {T,N} = :(result = Array{result_element_type}(undef, $(size(func_array))))

return_expression(::SArray) = :(return SArray(result))
return_expression(::Array) = :(return result)

_value(a::Node) = is_constant(a) ? value(a) : NaN
function _infer_numeric_eltype(array::AbstractArray{<:Node})
    eltype = Union{}
    for elt in array
        eltype = promote_type(eltype, typeof(_value(elt)))
    end
    eltype
end

"""Should only be called if `all_constants(func_array) == true`. Unpredictable results otherwise."""
function to_number(func_array::AbstractArray{T}) where {T<:Node}
    #find type
    element_type = _infer_numeric_eltype(func_array)
    tmp = similar(func_array, element_type)
    @. tmp = _value(func_array)
    return tmp
end

function to_number(func_array::SparseMatrixCSC{T}) where {T<:Node}
    nz = nonzeros(func_array)
    element_type = _infer_numeric_eltype(nz)
    tmp = similar(nz, element_type)
    @. tmp = value(nz)
    return tmp
end


"""
    make_Expr(
        func_array::AbstractArray{<:Node},
        input_variables::AbstractVector{<:Node},
        in_place::Bool,
        init_with_zeros::Bool
    )
"""
function make_Expr(func_array::AbstractArray{T}, input_variables::AbstractVector...; in_place::Bool=false, init_with_zeros::Bool=true) where {T<:Node}
    node_to_var = IdDict{Node,Union{Symbol,Real,Expr}}()
    body = Expr(:block)

    input_variable_names = Symbol[]

    node_to_index = IdDict{Node,Union{Expr,Int64}}()
    for (j, input_var_array) in pairs(input_variables)
        var_name = Symbol("input_variables$j")
        push!(input_variable_names, var_name)
        for (i, node) in pairs(input_var_array)
            node_to_index[node] = :($var_name[$i])
        end
    end

    num_zeros = count(is_identically_zero, (func_array))
    num_const = count((x) -> is_constant(x) && !is_identically_zero(x), func_array)


    zero_threshold = 0.5
    const_threshold = 0.5

    is_all_constant = num_const + num_zeros == length(func_array)
    is_mostly_zero = num_zeros > 5 * num_const

    # figure out if we have clear majority of terms and select the initialization strategy accordingly
    if (is_all_constant && is_mostly_zero) || (!is_all_constant && num_zeros > zero_threshold * length(func_array))
        initialization_strategy = :zero
    elseif (is_all_constant && !is_mostly_zero) || (!is_all_constant && num_const > const_threshold * length(func_array))
        initialization_strategy = :const
    else
        initialization_strategy = :undef
    end

    # declare result element type, and result variable if not provided by the user
    if in_place
        push!(body.args, :(result_element_type = promote_type(eltype.(($(input_variable_names...),))...)))
    else
        push!(body.args, :(result_element_type = promote_type($(_infer_numeric_eltype(func_array)), (eltype.(($(input_variable_names...),)))...)))
        push!(body.args, undef_array_declaration(func_array))
    end

    # if there was a clear majority of terms, initialize those in one shot to reduce code size
    if initialization_strategy === :zero && (init_with_zeros && in_place || !in_place)
        push!(body.args, :(result .= zero(result_element_type)))
    elseif initialization_strategy === :const
        push!(body.args, :(result .= $(to_number(func_array))))
    end



    for (i, node) in pairs(func_array)
        # skip all terms that we have computed above during construction
        if is_constant(node) && initialization_strategy === :const || # already initialized as constant above
           is_identically_zero(node) && (initialization_strategy === :zero || !init_with_zeros) # was already initialized as zero above or we don't want to initialize with zeros
            continue
        end
        node_body, variable = function_body!(node, node_to_index, node_to_var)
        for arg in node_body.args
            push!(body.args, arg)
        end
        push!(body.args, :(result[$i] = $variable))
    end

    # return result or nothing if in_place
    if in_place
        push!(body.args, :(return nothing))
    else
        push!(body.args, return_expression(func_array))
    end

    expected_input_lengths = map(length, input_variables)
    #boundscheck code

    input_check = :(
        @boundscheck begin

        lengths = $expected_input_lengths
        for (input, expected_length) in zip(($((input_variable_names)...),), lengths)
            if length(input) != expected_length
                actual_lengths = map(length, ($(input_variable_names...),))
                throw(ArgumentError("The input variables must have the same length as the input_variables argument to make_function. Expected lengths: $lengths. Actual lengths: $(actual_lengths)."))
            end
        end
    end)

    # wrap in function body
    if in_place
        expected_result_length = length(func_array)

        return :((result, $(input_variable_names...)) ->
            begin
                @boundscheck begin
                    expected_res_length = $expected_result_length
                    if length(result) != expected_res_length
                        throw(ArgumentError("The in place result vector does not have the expected length. Expected length: $expected_res_length. Actual length: $(length(result))."))
                    end
                end

                $input_check

                @inbounds begin
                    $body
                end
            end)
    else
        return :(($(input_variable_names...),) ->
            begin
                #here
                $input_check
                @inbounds begin
                    $body
                end
            end)
    end
end
export make_Expr

"""
    make_Expr(
        A::SparseMatrixCSC{<:Node,<:Integer},
        input_variables::AbstractVector{<:Node},
        in_place::Bool, init_with_zeros::Bool
    )

`init_with_zeros` argument is not used for sparse matrices."""
function make_Expr(A::SparseMatrixCSC{T,Ti}, input_variables::AbstractVector...; in_place::Bool=false, init_with_zeros::Bool=true) where {T<:Node,Ti}
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

    num_consts = count(x -> is_constant(x), vals)
    if num_consts == nnz(A) #all elements are constant
        push!(body.args, :(vals .= $(to_number(A))))
        if in_place
            return :((result, input_variables) -> $body)
        else
            push!(body.args, :(return result))
            return :((input_variables) -> $body)
        end
    else
        node_to_index = IdDict{Node,Union{Expr,Int64}}()
        for (i, node) in pairs(input_variables)
            node_to_index[node] = :(input_variables[$i])
        end

        for j = 1:n
            for i in nzrange(A, j)
                node_body, variable = function_body!(vals[i], node_to_index, node_to_var)
                for arg in node_body.args
                    push!(body.args, arg)
                end
                push!(node_body.args,)

                push!(body.args, :(vals[$i] = $variable))
            end
        end

        push!(body.args, :(return result))

        if in_place
            return :((result, input_variables::AbstractArray) -> $body)
        else
            return :((input_variables::AbstractArray) -> $body)
        end
    end
end
export make_Expr

"""
    make_function(
        func_array::AbstractArray{<:Node},
        input_variables::AbstractVector{<:Node}...;
        in_place::Bool=false, init_with_zeros::Bool=true
    )

Makes a function to evaluate the symbolic expressions in `func_array`. Every variable that is used in `func_array` must also be in `input_variables`. However, it will not cause an error if variables in `input_variables` are not variables used by `func_array`.

```julia
julia> @variables x
x

julia> f = x+1
(x + 1)


julia> jac = jacobian([f],[x]) #the Jacobian has a single constant element, 1, and is no longer a function of x
1×1 Matrix{FastDifferentiation.Node}:
 1

 julia> fjac = make_function(jac,[x])
 ...
 
 julia> fjac(2.0) #the value 2.0 is passed in for the variable x but has no effect on the output. Does not cause a runtime exception.
 1×1 Matrix{Float64}:
  1.0
```

If `in_place=false` then a new array will be created to hold the result each time the function is called. If `in_place=true` the function expects a user supplied array to hold the result. The user supplied array must be the first argument to the function.

```julia
julia> @variables x
x

julia> f! = make_function([x,x^2],[x],in_place=true)
...

julia> result = zeros(2)
2-element Vector{Float64}:
 0.0
 0.0

julia> f!(result,[2.0])
4.0

julia> result
2-element Vector{Float64}:
 2.0
 4.0
```

If the array is sparse then the keyword argument `init_with_zeros` has no effect. If the array is dense and `in_place=true` then the keyword argument `init_with_zeros` affects how the in place array is initialized. If `init_with_zeros = true` then the in place array is initialized with zeros. If `init_with_zeros=false` it is the user's responsibility to initialize the array with zeros before passing it to the runtime generated function.

This can be useful for modestly sparse dense matrices with say at least 1/4 of the array entries non-zero. In this case a sparse matrix may not be as efficient as a dense matrix. But a large fraction of time could be spent unnecessarily setting elements to zero. In this case you can initialize the in place Jacobian array once with zeros before calling the run time generated function.
"""
function make_function(func_array::AbstractArray{T}, input_variables::AbstractVector{<:Node}...; in_place::Bool=false, init_with_zeros::Bool=true) where {T<:Node}
    vars = variables(func_array) #all unique variables in func_array
    all_input_vars = vcat(input_variables...)

    #Because FD defines == operator for Node, which does not return a boolean, many builtin Julia functions will not work as expected. For example:
    #  vars ⊆ all_input_vars errors because internally issubset tests for equality between the node values using ==, not ===. == returns a Node value but the issubset function expects a Bool.

    temp = Vector{eltype(vars)}(undef, 0)

    input_dict = IdDict(zip(all_input_vars, all_input_vars))
    for one_var in vars
        value = get(input_dict, one_var, nothing)
        if value === nothing
            push!(temp, one_var)
        end
    end

    @assert length(temp) == 0 "The variables $temp were not in the input_variables argument to make_function. Every variable that is used in your function must have a corresponding entry in the input_variables argument."

    temp = make_Expr(func_array, input_variables..., in_place=in_place, init_with_zeros=init_with_zeros)
    @RuntimeGeneratedFunction(temp)
end
export make_function


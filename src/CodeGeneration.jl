function make_Expr(func_array::AbstractArray{T,N}, input_variables::AbstractVector{S}; in_place=false) where {T<:Node,S<:Node,N}
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

    if in_place
        return :((input_variables, result) -> $body)
    else
        return :((input_variables) -> $body)
    end
end
export make_Expr

make_function(func_array::AbstractArray{T}, input_variables::AbstractVector{S}; in_place=false, use_vector_runtime_args=false) where {T<:Node,S<:Node} = @RuntimeGeneratedFunction(make_Expr(func_array, input_variables, in_place=in_place))
export make_function
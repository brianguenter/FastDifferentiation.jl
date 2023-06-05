return_declaration(::StaticArray{S,T,N}, input_variables::AbstractVector{T2}) where {S,T,T2,N} = :(result = MArray{$(S),promote_type(Float64, eltype(input_variables)),$N}(undef))
return_declaration(func_array::Array{T,N}, input_variables::AbstractVector{S}) where {S,T,N} = :(result = Array{promote_type(Float64, eltype(input_variables)),$N}(undef, $(size(func_array)...)))

return_expression(::SArray) = :(return SArray(result))
return_expression(::Array) = :(return result)

function make_Expr(func_array::AbstractArray{T}, input_variables::AbstractVector{S}, in_place::Bool) where {T<:Node,S<:Node}
    node_to_var = Dict{Node,Union{Symbol,Real,Expr}}()
    body = Expr(:block)

    if !in_place
        push!(body.args, (return_declaration(func_array, input_variables)))
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

    if !in_place
        push!(body.args, return_expression(func_array))
    end

    if in_place
        return :((input_variables, result) -> $body)
    else
        return :((input_variables) -> $body)
    end
end
export make_Expr


# function update_sparse(mat::SparseMatrixCSC, partial_variables)
#need map from var_index index of function variable that is in this column
# for i in 1:length(mat.colptr)
#     f_index = mat.rowval[i]
#     v_index = mat.colptr[i]
#     nzval[i] = evaluate_path(gr, f_index, var_index) #not correct but closest
# end
function make_Expr(A::SparseMatrixCSC{T,Ti}, input_variables::AbstractVector{S}, in_place::Bool) where {T<:Node,S<:Node,Ti}
    rows = rowvals(A)
    vals = nonzeros(A)
    _, n = size(A)
    body = Expr(:block)
    node_to_var = Dict{Node,Union{Symbol,Real,Expr}}()

    if !in_place #have to store the sparse vector indices in the generated code to know how to create sparsity pattern
        push!(body.args, :(element_type = promote_type(Float64, eltype(input_variables))))
        push!(body.args, :(result = SparseMatrixCSC($(A.m), $(A.n), $(A.colptr), $(A.rowval), zeros(element_type, $(length(A.nzval))))))
    end

    push!(body.args, :(vals = nonzeros(result)))

    node_to_index = Dict{Node,Int64}()
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

make_function(func_array::SparseMatrixCSC{T,Ti}, input_variables::AbstractVector{S}; in_place=false) where {T<:Node,S<:Node,Ti} = @RuntimeGeneratedFunction(make_Expr(func_array, input_variables, in_place))

make_function(func_array::AbstractArray{T}, input_variables::AbstractVector{S}; in_place=false) where {T<:Node,S<:Node} = @RuntimeGeneratedFunction(make_Expr(func_array, input_variables, in_place))
export make_function
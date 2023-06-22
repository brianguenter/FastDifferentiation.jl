# struct FunctionCache{A<:AbstractArray{<:Node},N}
#     function_cache::Dict{BitVector,RuntimeGeneratedFunction}
#     original_function::A
#     inputs::NTuple{Vector{Node}}
#     operation::Function #or could do Union of the various functions that are legal
#     function_arguments::B #vector(vector(node)) but type may change
#     in_place::Bool

#     FunctionCache(orig_func::AbstractArray{<:Node}, operation::Function, function_args::AbstractVector{<:AbstractVector{<:Node}}, in_place::Bool)
# end

# (a::FunctionCache)(inputs::AbstractVector{T}) where {T<:Real}
# function make_function(func_array::AbstractArray{<:Node}, inputs::AbstractVector{<:Node}..., in_place=true)
#     if has_ifelse(func_array)
#         _make_conditional_function(func_array, inputs, in_place)
#     else
#         _make_function(func_array, inputs, in_place) #this is what I currently have.
#     end
# end

# function _make_conditional_function(func_array::T, inputs::AbstractVector{<:Node}..., in_place=true) where {T<:AbstractArray{<:Node}}
#     cache = FunctionCache{BitVector,T}()
#     return FunctionCache{T,length(inputs)}(cache, func_array, inputs)

# end

# is_conditional(a::Node) = is_tree(a) && value(a) in BOOL_OPS

# ifelse_nodes(nodes::Vector{Node}) = filter(x -> x === ifelse, nodes) #nodes(gr) is sorted by postorder number for results are as well.

# bool_nodes(nodes::Vector{Node}) = filter(x -> x in BOOL_OPS, nodes) #nodes(gr) is sorted by postorder number for results are as well. 



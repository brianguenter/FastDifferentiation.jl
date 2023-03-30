function_array(a::Symbol, n::Int, var::T) where {T<:Node} = [FastSymbolicDifferentiation.function_of(Symbol(a, x), var) for x in 1:n]
export function_array

random_joints(joint_functions::AbstractVector{<:Node}) = [transformation(random_rotation(x), SVector{3}(Node.(rand(3)))) for x in joint_functions]
export random_joints

struct Joints
    Ai::Vector{Matrix{Node}}

    Joints(joint_functions::AbstractVector{<:Node}) = new(random_joints(joint_functions))
end
export Joints

W(linkage::Joints, index::Int) = reduce(*, linkage.Ai[1:index], init=I())
export W

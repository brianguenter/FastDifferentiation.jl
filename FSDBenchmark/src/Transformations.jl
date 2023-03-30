# function quaternion(x::T, y::T, z::T) where {T<:Real}
#     return SVector{4,T}(sqrt(1 - (wx * wx - wy * wy - wz * wz / 4)), wx / 2, wy / 2, wz / 2)
# end

# function quaternion_multiplication(q::SVector{T,4}, p::SVector{S,4}) where {T<:Real,S<:Real}
#     a, b, c, d = q
#     e, f, g, h = p
#     return SVector{4,promote_type(T, S)}(a * e - b * f - c * g - d * h, a * f + b * e + c * h - d * g, a * g - b * h + c * e + d * f, a * h + b * g - c * f + d * e)
# end

# conjugate(a::SVector{T,4}) where {T} = SVector{T,4}(a[1], -a[2], -a[3], -a[4])
# inverse(a::SVector{T,4}) where {T} = conjugate(a) / (a ⋅ a) #only need to normalize if not unit quaternions

# function quaternion_rotation(q::SVector{T,4}, p::SVector{S,3}) where {T<:Real,S<:Real}
#     p′ = SVector{S,4}(p[1], p[2], p[3], zero(S))
#     return quaternion_multiplication(q, quaternion_multiplication(p′, conjugate(q)))
# end

import Base

cross(a::AbstractVector{T}) where {T} = SMatrix{3,3,T}(0, a[4], -a[3], -a[4], 0, a[2], a[3], -a[2], 0)
export cross

function normalize(a::AbstractVector)
    tmp = mapreduce(x -> x^2, +, a)
    return map(x -> x / tmp, a)
end

angle_axis(a::AbstractVector{T}) where {T} = SVector{4,T}(a[1], normalize(a[2:4])...)
export angle_axis

random_rotation(a::Node{<:FastSymbolicDifferentiation.UnspecifiedFunction}) = angle_axis([a, rand(3)...])
export random_rotation

"""convert angle axis to matrix form"""
function matrix(a::AbstractVector)
    @assert length(a) == 4
    return [1 0 0; 0 1 0; 0 0 1] + map(x -> x * sin(a[1]), cross(a)) + map(x -> (1 - cos(a[1])) * x, (cross(a) * cross(a)))
end
export matrix

function transformation(rotation::SVector{4,<:Node}, translation::SVector{3,<:Node})
    result = Matrix{Node}(undef, 4, 4)
    copyto!(result, CartesianIndices((1:3, 1:3)), matrix(angle_axis(rotation)), CartesianIndices((1:3, 1:3)))

    result[1, 4] = translation[1]
    result[2, 4] = translation[2]
    result[3, 4] = translation[3]
    result[4, 1] = 0.0
    result[4, 2] = 0.0
    result[4, 3] = 0.0
    result[4, 4] = 1

    return SMatrix{4,4}(result)
end
export transformation

I() = SMatrix{4,4}(
    1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0
)


using LinearAlgebra

function quaternion(x::T, y::T, z::T) where {T<:Real}
    return SVector{4,T}(sqrt(1 - (wx * wx - wy * wy - wz * wz / 4)), wx / 2, wy / 2, wz / 2)
end

function quaternion_multiplication(q::SVector{T,4}, p::SVector{S,4}) where {T<:Real,S<:Real}
    a, b, c, d = q
    e, f, g, h = p
    return SVector{4,promote_type(T, S)}(a * e - b * f - c * g - d * h, a * f + b * e + c * h - d * g, a * g - b * h + c * e + d * f, a * h + b * g - c * f + d * e)
end

conjugate(a::SVector{T,4}) where {T} = SVector{T,4}(a[1], -a[2], -a[3], -a[4])
inverse(a::SVector{T,4}) where {T} = conjugate(a) / (a ⋅ a) #only need to normalize if not unit quaternions

function quaternion_rotation(q::SVector{T,4}, p::SVector{S,3}) where {T<:Real,S<:Real}
    p′ = SVector{S,4}(p[1], p[2], p[3], zero(S))
    return quaternion_multiplication(q, quaternion_multiplication(p′, conjugate(q)))
end

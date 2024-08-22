function all_combinations(n::Integer)
    num = 2^n

    result = Vector{Vector{Bool}}(undef, num)

    for i in 0:num-1
        result[i+1] = Bool.(digits(i, base=2, pad=n))
    end
    return result
end
export all_combinations

struct Combinations{T<:Integer}
    n::T
end

iterate(a::Combinations) = (a, 0)

function iterate(a::Combinations{T}, state::T) where {T<:Integer}
    if state == 2^a.n
        return nothing
    else
        return (BitVector(Bool.(digits(state, base=2, pad=a.n))), state + 1)
    end
end

# Base.getindex(a::Combinations{T}, b::T) where {T<:Integer} = b
# length(a::Combinations) = 2^a.n
# eltype(a::Combinations) = Combinations
# size(a::Combinations, dims...) = length(a)
#These functions appear to be """returns true if a ⊆ b, false otherwise"""
bitvector_cache::Vector{BitVector} = BitVector[]

"""WARNING not multithread safe"""
function subset(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    index = findfirst(x -> length(x) == length(a), bitvector_cache)
    if index !== nothing
        temp = bitvector_cache[index]
    else
        temp = BitVector(undef, length(a))
        push!(bitvector_cache, temp)
    end
    temp .= a .& b
    sum(temp) == sum(a)
end
export subset

"""returns the set that has all the elements of b removed from a. This allocates"""
function set_diff(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    return @. !(a & b) & a
end
export set_diff

function bit_equal(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    return !any(a .⊻ b)
end
export bit_equal

is_zero(a::BitVector) = !any(a)
export is_zero

function to_string(a::BitVector, b::String)
    doms = ""
    for (i, bitval) in pairs(a)
        if bitval
            doms *= "$b$i,"
        end
    end
    return doms
end
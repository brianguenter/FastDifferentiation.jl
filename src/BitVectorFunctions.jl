#These functions appear to be 
bitvector_cache::Vector{BitVector} = BitVector[]

function get_bitvector(numbits::T) where {T<:Integer}
    global bitvector_cache

    index = findfirst(x -> length(x) == numbits, bitvector_cache)
    if index !== nothing
        temp = bitvector_cache[index]
    else
        temp = BitVector(undef, length(a))
        push!(bitvector_cache, temp)
    end
    return temp
end

function reclaim_bitvector(a::BitVector)
    global bitvector_cache

    push!(bitvector_cache, a)
end

"""WARNING not multithread safe

returns true if a ⊆ b, false otherwise"""
function subset(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    temp = get_bitvector(length(a))
    temp .= a .& b
    result = sum(temp) == sum(a)
    reclaim_bitvector(temp)
    return result
end
export subset

"""returns the set that has all the elements of b removed from a. This allocates"""
function set_diff(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    return @. !(a & b) & a
end
export set_diff

function overlap(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    temp = get_bitvector(length(a))
    temp .= a .& b
    result = any(temp)
    reclaim_bitvector(temp)
    return result
end

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
#These functions appear to be 
bitvector_cache::Vector{BitVector} = BitVector[]

peak_bitvector_cache_size::Int64 = 0

function get_bitvector(numbits::T) where {T<:Integer}
    global bitvector_cache
    global peak_bitvector_cache_size

    index = findfirst(x -> length(x) == numbits, bitvector_cache)
    if index !== nothing
        temp = bitvector_cache[index]
    else
        temp = BitVector(undef, numbits)
        push!(bitvector_cache, temp)
        if length(bitvector_cache) > peak_bitvector_cache_size
            peak_bitvector_cache_size = length(bitvector_cache)
        end
    end
    return temp
end


"""WARNING not multithread safe

returns true if a ⊆ b, false otherwise"""
function subset(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    temp = get_bitvector(length(a))
    temp .= a .& b
    result = sum(temp) == sum(a)
    return result
end


"""returns the set that has all the elements of b removed from a. This allocates"""
function set_diff(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    return @. !(a & b) & a
end


"""removes elements of b from a"""
function set_diff!(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    @. a = !(a & b) & a
    return nothing
end


"""returns true if a ∩ b is not empty"""
function overlap(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    temp = get_bitvector(length(a))
    temp .= a .& b
    result = any(temp)
    return result
end

function bit_equal(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    return !any(a .⊻ b)
end

is_zero(a::BitVector) = !any(a)


function to_string(a::BitVector, b::String)
    doms = ""
    for (i, bitval) in pairs(a)
        if bitval
            doms *= "$b$i,"
        end
    end
    return doms
end
#These functions appear to be 
bitvector_cache::Dict{Int64,BitVector} = Dict{Int64,BitVector}()

peak_bitvector_cache_size::Int64 = 0

#Really should occasionally 

function get_bitvector(numbits::T) where {T<:Integer}
    global bitvector_cache
    global peak_bitvector_cache_size

    tmp = get(bitvector_cache, numbits, nothing)
    if tmp !== nothing
        temp = bitvector_cache[numbits]
    else
        temp = BitVector(undef, numbits)
        bitvector_cache[numbits] = temp
    end
    return temp
end


"""
    subset(a::BitVector, b::BitVector)

Returns true if a ⊆ b, false otherwise.
WARNING not multithread safe."""
function subset(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    temp = get_bitvector(length(a))
    temp .= a .& b
    result = sum(temp) == sum(a)
    return result
end


"""
    set_diff(a::BitVector, b::BitVector)

Returns the set that has all the elements of b removed from a.
This allocates."""
function set_diff(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    return @. !(a & b) & a
end


"""
    set_diff!(a::BitVector, b::BitVector)

Removes elements of b from a."""
function set_diff!(a::BitVector, b::BitVector)
    @assert length(a) == length(b)
    @. a = !(a & b) & a
    return nothing
end


"""
    overlap(a::BitVector, b::BitVector)

Returns true if a ∩ b is not empty."""
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

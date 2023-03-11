"""returns true if a ⊆ b, false otherwise"""
subset(a::BitVector, b::BitVector) = sum(a .& b) == sum(a)
export subset

"""returns the set that has all the elements of b removed from a"""
set_diff(a::BitVector, b::BitVector) = @. !(a & b) & a
export set_diff

set_zero!(a::BitVector) = @. a = a ⊻ a
export zero_out

bit_equal(a::BitVector, b::BitVector) = !any(a .⊻ b)
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
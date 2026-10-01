function all_combinations(n::Integer)
    num = 2^n

    result = Vector{Vector{Bool}}(undef, num)

    for i in 0:num-1
        result[i+1] = Bool.(digits(i, base=2, pad=n))
    end
    return result
end
export all_combinations

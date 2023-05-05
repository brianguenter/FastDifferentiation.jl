@memoize function Chebyshev(n, x)
    if n == 0
        return 1
    elseif n == 1
        return x
    else
        return 2 * (x) * Chebyshev(n - 1, x) - Chebyshev(n - 2, x)
    end
end
export Chebyshev

@memoize function Chebyshev(n, x::Node{T,0}) where {T}
    if n == 0
        return 1
    elseif n == 1
        return x
    else
        return 2 * (x) * Chebyshev(n - 1, x) - Chebyshev(n - 2, x)
    end
end
export Chebyshev

"""WARNING: slow for large n"""
Chebyshev_exe(n, x::Symbolics.Num) = make_function(DerivativeGraph([expr_to_dag(Chebyshev(n, x))]))
export Chebyshev_exe

Chebyshev_exe(n, x::Node) = make_function(DerivativeGraph([Chebyshev(n, x)]))

function chebyshev(model_size)
    Symbolics.@variables x
    nx = Node(x)
    # return RnToRmGraph([expr_to_dag(Chebyshev(order, x))])
    return DerivativeGraph([Chebyshev(model_size, nx)])
end
export chebyshev


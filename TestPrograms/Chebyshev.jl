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

@memoize function Chebyshev(n, x::Node)
    @assert is_variable(x)

    if n == 0
        return 1
    elseif n == 1
        return x
    else
        return 2 * (x) * Chebyshev(n - 1, x) - Chebyshev(n - 2, x)
    end
end
export Chebyshev

Chebyshev_exe(n, x::Node) = make_function(DerivativeGraph([Chebyshev(n, x)]))


function chebyshev(::FastSymbolic, model_size)
    @variables x

    # return RnToRmGraph([expr_to_dag(Chebyshev(order, x))])
    return DerivativeGraph([Chebyshev(model_size, x)])
end
export chebyshev


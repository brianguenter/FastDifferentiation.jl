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

@memoize function Chebyshev(n, x::FD.Node)
    @assert FD.is_variable(x)

    if n == 0
        return 1
    elseif n == 1
        return x
    else
        return 2 * (x) * Chebyshev(n - 1, x) - Chebyshev(n - 2, x)
    end
end
export Chebyshev

Chebyshev_exe(n, x::FD.Node) = make_function(FD.DerivativeGraph([Chebyshev(n, x)]))


function chebyshev(::FastSymbolic, model_size)
    @variables x

    # return RnToRmGraph([expr_to_dag(Chebyshev(order, x))])
    return FD.DerivativeGraph([Chebyshev(model_size, x)])
end
export chebyshev


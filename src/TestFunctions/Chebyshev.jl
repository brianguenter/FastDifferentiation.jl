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
Chebyshev_exe(n, x::Symbolics.Num) = to_function(RnToRmGraph([expr_to_dag(Chebyshev(n, x))]))
export Chebyshev_exe

Chebyshev_exe(n, x::Node) = to_function(RnToRmGraph([Chebyshev(n, x)]))

function chebyshev_graph(order)
    Symbolics.@variables x
    nx = Node(x)
    # return RnToRmGraph([expr_to_dag(Chebyshev(order, x))])
    return RnToRmGraph([Chebyshev(order, nx)])
end
export chebyshev_graph


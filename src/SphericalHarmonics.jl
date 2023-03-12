module SphericalHarmonics
using Symbolics
using FastSymbolicDifferentiation
using Memoize

@memoize function P(l, m, z)
    if l == 0 && m == 0
        return 1.0
    elseif l == m
        return (1 - 2m) * P(m - 1, m - 1, z)
    elseif l == m + 1
        return (2m + 1) * z * P(m, m, z)
    else
        return ((2l - 1) / (l - m) * z * P(l - 1, m, z) - (l + m - 1) / (l - m) * P(l - 2, m, z))
    end
end
export P

@memoize function S(m, x, y)
    if m == 0
        return 0
    else
        return x * C(m - 1, x, y) - y * S(m - 1, x, y)
    end
end
export S

@memoize function C(m, x, y)
    if m == 0
        return 1
    else
        return x * S(m - 1, x, y) + y * C(m - 1, x, y)
    end
end
export C

@memoize function N(l, m)
    @assert m >= 0
    if m == 0
        return sqrt((2l + 1 / (4π)))
    else
        # return sqrt((2l+1)/2π * factorial(big(l-m))/factorial(big(l+m)))
        #use factorial_approximation instead of factorial because the latter does not use Stirlings approximation for large n. Get error for n > 2 unless using BigInt but if use BigInt get lots of rational numbers in symbolic result.
        return sqrt((2l + 1) / 2π * factorial_approximation(l - m) / factorial_approximation(l + m))
    end
end
export N

"""l is the order of the spherical harmonic. I think"""
@memoize function Y(l, m, x, y, z)
    @assert l >= 0
    @assert abs(m) <= l
    if m < 0
        return N(l, abs(m)) * P(l, abs(m), z) * S(abs(m), x, y)
    else
        return N(l, m) * P(l, m, z) * C(m, x, y)
    end
end
export Y

SHFunctions(max_l, x::Node, y::Node, z::Node) = SHFunctions(Vector{Node}(undef, 0), max_l, x, y, z)
SHFunctions(max_l, x::Num, y::Num, z::Num) = SHFunctions(Vector{Num}(undef, 0), max_l, x, y, z)

function SHFunctions(shfunc, max_l, x, y, z)
    for l in 0:max_l-1
        for m in -l:l
            push!(shfunc, Y(l, m, x, y, z))
        end
    end

    return shfunc
end
export SHFunctions


function SHDerivatives(max_l, x::Num, y::Num, z::Num)
    local derivatives = Vector{Num}(undef, 0)
    dx = Differential(x)
    dy = Differential(y)
    dz = Differential(z)

    for l in 0:max_l-1
        for m in -l:l
            tmp = Y(l, m, x, y, z)
            push!(derivatives, expand_derivatives(dx(tmp), true))
            push!(derivatives, expand_derivatives(dy(tmp), true))
            push!(derivatives, expand_derivatives(dz(tmp), true))
        end
    end

    return derivatives
end
export SHDerivatives

to_dag(max_l, x, y, z) = expr_to_dag.(SHDerivatives(max_l, x, y, z))
export to_dag

function to_graph(max_l)
    @variables x, y, z
    nx = Node(x)
    ny = Node(y)
    nz = Node(z)

    graph = RnToRmGraph(SHFunctions(max_l, nx, ny, nz))
    return graph, x, y, z
end
export to_graph

function spherical_harmonics_jacobian(order)
    gr, x, y, z = to_graph(order)
    return jacobian_function!(gr)
end
export spherical_harmonics_jacobian

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

Chebyshev_exe(n, x) = to_function(RnToRmGraph([expr_to_dag(Chebyshev(n, x))]))
export Chebyshev_exe

function chebyshev_graph(order)
    @variables x
    nx = Node(x)
    # return RnToRmGraph([expr_to_dag(Chebyshev(order, x))])
    return RnToRmGraph([Chebyshev(order, nx)])
end
export chebyshev_graph

end #module
export SphericalHarmonics

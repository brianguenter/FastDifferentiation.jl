# Symbolic Processing

Because **FD** can generate true symbolic derivatives it can easily be used in conjunction with Symbolics.jl using the package [FDConversion.jl](https://github.com/brianguenter/FDConversion.jl) (still under development).

A rule of thumb is that if your function is small (a few hundred operations or less) or tree like (where each node in the expression graph has one parent on average) then Symbolics.jl may outperform or equal **FD**. For more complex functions with many common subexpressions **FD** may substantially outperform Symbolics.jl.
 
Take these benchmarks with a large grain of salt since there are so few of them. Whether your function will have this kind of performance improvement relative to Symbolics.jl is hard to predict until the benchmark set gets much bigger.

These benchmarks should give you a sense of what performance you might achieve for symbolic processing. There are three types of benchmarks: **Symbolic**, **MakeFunction**, and **Exe**.

* The **Symbolic** benchmark is the time required to compute just the symbolic form of the derivative. The Symbolic benchmark can be run with simplification turned on or off for Symbolics.jl. If simplification is on then computation time can be extremely long but the resulting expression might be simpler and faster to execute.

* The **MakeFunction** benchmark is the time to generate a Julia Expr from an already computed symbolic derivative and to then compile it.

* The **Exe** benchmark measures just the time required to execute the compiled function using an in-place matrix.

All benchmarks show the ratio of time taken by Symbolics.jl to FastDifferentiation.jl. Numbers greater than 1 mean FastDifferentiation is faster.

All benchmarks were run on an AMD Ryzen 9 7950X 16-Core Processor with 32GB RAM running Windows 11 OS, Julia version 1.9.0.
### Chebyshev polynomial
The first example is a recursive function for 
the Chebyshev polynomial of order n:

```
@memoize function Chebyshev(n, x)
    if n == 0
        return 1
    elseif n == 1
        return x
    else
        return 2 * (x) * Chebyshev(n - 1, x) - Chebyshev(n - 2, x)
    end
end
```
The function is memoized so the recursion executes efficiently. 

The recursive function returns an nth order polynomial in the variable x. The derivative of this polynomial would be order n-1 so a perfect symbolic simplification would result in a function with 2*(n-2) operations. For small values of n Symbolics.jl simplification does fairly well but larger values result in very inefficient expressions.

Because **FD** doesn't do sophisticated symbolic simplification it generates a derivative with approximately 2.4x the number of operations in the original recursive expression regardless of n. This is a case where a good hand generated derivative would be more efficient than **FD**.

The Chebyshev expression graph does not have many nodes even at the largest size tested (graph size increases linearly with Chebyshev order).

The first set of three benchmarks show results with simplification turned off in Symbolics.jl, followed by a set of three with simplification turned on. Performance is somewhat better in the latter case but still slower than the FD executable. Note that the y axis is logarithmic.

#### Chebyshev benchmarks with simplification off
![Symbolic processing, simplify=false](Illustrations/figure_chebyshev_Symbolic_simplify_false.svg) 
![MakeFunction, simplify=false](Illustrations/figure_chebyshev_MakeFunction_simplify_false.svg) 
![Exe, simplify=false](Illustrations/figure_chebyshev_Exe_simplify_false.svg)



#### Chebyshev benchmarks with simplification on
![MakeFunction, simplify=false](Illustrations/figure_chebyshev_Exe_simplify_true.svg)

With simplification on performance of the executable derivative function for Symbolics.jl is slightly better than with simplification off. But simplification processing time is longer.
 
### Spherical Harmonics

The second example is the spherical harmonics function. This is the expression graph for the spherical harmonic function of order 8:
![MakeFunction, simplify=false](Illustrations/sphericalharmonics_8.svg)

```
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

function factorial_approximation(x)
    local n1 = x
    sqrt(2 * π * n1) * (n1 / ℯ * sqrt(n1 * sinh(1 / n1) + 1 / (810 * n1^6)))^n1
end
export factorial_approximation

function compare_factorial_approximation()
    for n in 1:30
        println("n $n relative error $((factorial(big(n))-factorial_approximation(n))/factorial(big(n)))")
    end
end
export compare_factorial_approximation

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

"""l is the order of the spherical harmonic"""
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
SHFunctions(max_l, x::Symbolics.Num, y::Symbolics.Num, z::Symbolics.Num) = SHFunctions(Vector{Symbolics.Num}(undef, 0), max_l, x, y, z)

function SHFunctions(shfunc, max_l, x, y, z)
    for l in 0:max_l-1
        for m in -l:l
            push!(shfunc, Y(l, m, x, y, z))
        end
    end

    return shfunc
end
export SHFunctions

function spherical_harmonics(::JuliaSymbolics, model_size)
    Symbolics.@variables x y z
    return SHFunctions(model_size, x, y, z), [x, y, z]
end

function spherical_harmonics(::FastSymbolic, model_size, x, y, z)
    graph = DerivativeGraph(SHFunctions(model_size, x, y, z))
    return graph
end

function spherical_harmonics(package::FastSymbolic, model_size)
    FD.@variables x, y, z
    return spherical_harmonics(package, model_size, x, y, z)
end
export spherical_harmonics
```

As was the case for Chebyshev polynomials the number of paths from the roots to the variables is much greater than the number of nodes in the graph. Once again the y axis is logarithmic.

![Symbolic processing, simplify=false](Illustrations/figure_spherical_harmonics_Symbolic_simplify_false.svg)
![MakeFunction, simplify=false](Illustrations/figure_spherical_harmonics_MakeFunction_simplify_false.svg)
![Exe, simplify=false](Illustrations/figure_spherical_harmonics_Exe_simplify_false.svg)
 
 The **Exe** benchmark took many hours to run and was stopped at model size 24 instead of 25 as for the **Symbolic** and **MakeFunction** benchmarks.

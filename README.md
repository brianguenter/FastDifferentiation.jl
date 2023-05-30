# FastDifferentiation

[![Build Status](https://github.com/brianguenter/FastDifferentiation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brianguenter/FastDifferentiation.jl/actions/workflows/CI.yml?query=branch%3Amain)


FastDifferentiation (**FD**) is a package for generating efficient executables to evaluate derivatives of Julia functions. It can also generate efficient true symbolic derivatives for symbolic analysis. Unlike forward and reverse mode automatic differentiation **FD** automatically generates efficient derivatives for arbitrary function types: ℝ¹->ℝ¹, ℝ¹->ℝᵐ, ℝⁿ->ℝ¹, and ℝⁿ->ℝᵐ, m≠1,n≠1. 

Compared to forward and reverse mode AD the compiled **FD** derivative executable may have somewhat better performance for the ℝ¹->ℝ¹, ℝ¹->ℝᵐ, and ℝⁿ->ℝ¹ function types. This is because the runtime overhead present in most forward and reverse AD algorithms is compiled out in **FD**. 

For the ℝⁿ->ℝᵐ case the performance improvement relative to conventional AD algorithms may be larger. The **FD** algorithm finds expressions shared between partials and computes them only once - in this case **FD** executables may be significantly more efficient than either forward or reverse mode AD. In some cases **FD** derivatives can be as efficient as manually coded derivatives[^d].

 **FD** may take much less time to compute symbolic derivatives than Symbolics.jl even in the ℝ¹->ℝ¹ case. The executables generated by **FD** may also be much faster (see [Symbolic Processing](#SymbolicProcessing)). 

You should consider using FastDifferentiation when you need: 
* a fast executable for evaluating the derivative of a function and the overhead of the preprocessing/compilation time is swamped by evaluation time.
* additional symbolic processing on your derivative. **FD** can generate a true symbolic derivative to be processed further in Symbolics.jl or another computer algebra system (CAS).

This is the **FD** feature set (operations marked ❌ will be completed soon):

<table>
<tr>
<td> <b></b>
<td> <b>Dense Jacobian</b> <td>  <b>Sparse Jacobian</b> </td> 
<td>  <b>Dense Hessian</b> </td><td>  <b> Sparse Hessian</b> </td> 
<td>  <b>Higher order derivatives</b> </td> 
<td>  <b>Jᵀv</b> </td> 
<td>  <b>Jv</b> </td> 
<td> <b> Hv </b> </td>
</tr>
<tr>
<td> <b> Compiled function </b> </td> 
<td> ✅ </td>
<td> ❌ </td>
<td> ✅ </td>
<td> ❌  </td>
<td> ✅ </td>
<td> ✅ </td>
<td> ✅ </td>
<td> ✅ </td>
</tr>
<tr>
<td> <b> Symbolic expression </b> </td> 
<td> ✅ </td>
<td> ✅ </td>
<td> ✅ </td>
<td> ❌  </td>
<td> ✅ </td>
<td> ✅ </td>
<td> ✅ </td>
<td> ✅ </td>
</tr>

</table>

Jᵀv and Jv compute the Jacobian transpose times a vector and the Jacobian times a vector, without explicitly forming the Jacobian matrix. For applications see this [paper](https://arxiv.org/abs/1812.01892). 

Hv computes the Hessian times a vector without explicitly forming the Hessian matrix. This can be useful when the Hessian matrix is large and sparse.

If you use FD in your work please share the functions you differentiate with me. I'll add them to the benchmarks. The more functions available to test the easier it is for others to determine if FD will help with their problem.

This is **beta** software being modified on a daily basis. Expect bugs and frequent, possibly breaking changes, over the next month or so. Documentation is frequently updated so check the latest docs before filing an issue. Your problem may have been fixed and documented.

<details> 
 <summary> <b> Examples and basic usage </b> </summary>
 
The first step is to create **FD** variables which are then passed to the function you want to differentiate. The return value is a graph structure which **FD** will analyze to generate efficient executables or symbolic expressions.
 
**FD** uses a global cache for common subexpression elimination so the **FD** expression preprocessing step is not thread safe. 

Under ordinary conditions the memory used by the cache won't be an issue. But, if you have a long session where you are creating many complex functions it is possible the cache will use too much memory. If this happens call the function `clear_cache` after you have completely processed your expression.

Set up variables:
```
using FastDifferentiation

@variables x y z

```
 Make a vector of variables
 ```
julia> X = make_variables(:x,3)
3-element Vector{Node}:
 x1
 x2
 x3
```
Make an executable function
```
julia> xy_exe = make_function([x^2*y^2,sqrt(x*y)],[x,y]) #[x,y] vector specifies the order of the arguments to the exe
...

julia> xy_exe([1.0,2.0])
2-element Vector{Float64}:
 4.0
 1.4142135623730951
 ```
Compute Hessian:
```
@variables x y z

julia> h_symb = hessian(x^2+y^2+z^2,[x,y,z])
3×3 Matrix{Node}:
 2    0.0  0.0
 0.0  2    0.0
 0.0  0.0  2

julia> h_symb1 = hessian(x^2*y^2*z^2,[x,y,z])
3×3 Matrix{FastDifferentiation.Node}:
 (2 * ((z ^ 2) * (y ^ 2)))        (((2 * x) * (2 * y)) * (z ^ 2))  (((2 * x) * (2 * z)) * (y ^ 2))
 (((2 * y) * (2 * x)) * (z ^ 2))  (2 * ((z ^ 2) * (x ^ 2)))        (((2 * y) * (2 * z)) * (x ^ 2))
 (((2 * z) * (2 * x)) * (y ^ 2))  (((2 * z) * (2 * y)) * (x ^ 2))  (2 * ((x ^ 2) * (y ^ 2)))

julia> hexe_1 = make_function(h_symb1,[x,y,z])
...
julia> hexe_1([1.0,2.0,3.0])
3×3 Matrix{Float64}:
 72.0  72.0  48.0
 72.0  18.0  24.0
 48.0  24.0   8.0
```
Compute `Hv` without forming the full Hessian matrix. This is useful if the Hessian is very large
```
julia> @variables x y
y

julia> f = x^2 * y^2
((x ^ 2) * (y ^ 2))

julia> hv_fast, v_vec2 = hessian_times_v(f, [x, y])
...

julia> hv_fast_exe = make_function(hv_fast, [[x, y]; v_vec2]) #need v_vec2 because hv_fast is a function of x,y,v1,v2 and have to specify the order of all inputs to the executable
...
julia> hv_fast_exe([1.0,2.0,3.0,4.0]) #first two vector elements are x,y last two are v1,v2
2-element Vector{Float64}:
 56.0
 32.0
```
Compute Jacobian:
```
julia> f1 = cos(x) * y
(cos(x) * y)

julia> f2 = sin(y) * x
(sin(y) * x)

julia> symb = jacobian([f1, f2], [x, y]) #non-destructive
2×2 Matrix{Node}:
 (y * -(sin(x)))  cos(x)
 sin(y)           (x * cos(y))
```
Create executable to evaluate Jacobian:
```
julia> jac_exe = make_function(symb,[x,y])
...
julia> jac_exe([1.0,2.0])
2×2 Matrix{Float64}:
 -1.68294    0.540302
  0.909297  -0.416147
```
Executable with in_place matrix evaluation to avoid allocation of a matrix for the Jacobian (in_place option available on all executables including Jᵀv,Jv,Hv):
```
julia> jac_exe = make_function(symb,[x,y], in_place=true)
...
julia> a = Matrix{Float64}(undef,2,2)
2×2 Matrix{Float64}:
 0.0  0.0
 0.0  6.93532e-310

julia> jac_exe([1.0,2.0],a)
2×2 Matrix{Float64}:
 -1.68294    0.540302
  0.909297  -0.416147

julia> a
2×2 Matrix{Float64}:
 -1.68294    0.540302
  0.909297  -0.416147
```

For faster execution call the executable function with an `SVector` (for short vectors, probably < 100 elements):
```
julia> jac_exe(SVector{2}(1.0,2.0))
2×2 Matrix{Float64}:
 -1.68294    0.540302
  0.909297  -0.416147
 ```
Compute any subset of the columns of the Jacobian:
```
julia> symb = jacobian([x*y,y*z,x*z],[x,y,z]) #all columns
3×3 Matrix{Node}:
 y    x    0.0
 0.0  z    y
 z    0.0  x

julia> symb = jacobian([x*y,y*z,x*z],[x,y]) #first two columns
3×2 Matrix{Node}:
 y    x
 0.0  z
 z    0.0

julia> symb = jacobian([x*y,y*z,x*z],[z,y]) #second and third columns, reversed so ∂f/∂z is 1st column of the output, ∂f/∂y the 2nd
3×2 Matrix{Node}:
 0.0  x
 y    z
 x    0.0
 ```

Symbolic and executable Jᵀv and Jv (see this [paper](https://arxiv.org/abs/1812.01892) for applications of this operation).
```
julia> (f1,f2) = cos(x)*y,sin(y)*x
((cos(x) * y), (sin(y) * x))

julia> jv,vvec = jacobian_times_v([f1,f2],[x,y])
...

julia> jv_exe = make_function(jv,[[x,y];vvec])
...

julia> jv_exe([1.0,2.0,3.0,4.0]) #first 2 arguments are x,y values and last two are v vector values

2×1 Matrix{Float64}:
 -2.8876166853748195
  1.0633049342884753

julia> jTv,rvec = jacobian_transpose_v([f1,f2],[x,y])
...

julia> jtv_exe = make_function(jTv,[[x,y];rvec])
...
julia> jtv_exe([1.0,2.0,3.0,4.0])
2-element Vector{Float64}:
 -1.4116362015446517
 -0.04368042858415033
```

Convert between FastDifferentiation and Symbolics representations (requires [FDConversion](https://github.com/brianguenter/FDConversion/tree/main) package, not released yet[^b]):
```
julia> f = x^2+y^2 #Symbolics expression
x^2 + y^2

julia> Node(f) #convert to FastDifferentiation form
x^2 + y^2

julia> typeof(ans)
Node{SymbolicUtils.BasicSymbolic{Real}, 0}

julia> node_exp = x^3/y^4 #FastDifferentiation expression
((x ^ 3) / (y ^ 4))

julia> to_symbolics(node_exp)
(x^3) / (y^4)

julia> typeof(ans)
Symbolics.Num
```
</details>

<div id="SymbolicProcessing"></div>

<details>
    <summary> <b> Symbolic Processing </b> </summary>
 
##  Symbolic Processing

Because **FD** can generate true symbolic derivatives it can easily be used in conjunction with Symbolics.jl.

A rule of thumb is that if your function is small (a few hundred operations or less) or tree like (where each node in the expression graph has one parent on average) then Symbolics.jl may outperform or equal **FD**. For more complex functions with many common subexpressions **FD** may substantially outperform Symbolics.jl.
 
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

The Chebyshev expression graph does not have many nodes even at the largest size tested (graph size increases linearly with Chebyshev order). For example, here is the graph of the 10th order expression: 
<img src="Illustrations/chebyshev10.svg" alt="drawing" height="400">
The complexity arises from the number of different paths from the root to the leaf of the graph.

The first set of three benchmarks show results with simplification turned off in Symbolics.jl, followed by a set of three with simplification turned on. Performance is somewhat better in the latter case but still slower than the FD executable. Note that the y axis is logarithmic.

#### Chebyshev benchmarks with simplification off
<img src="Illustrations\figure_chebyshev_Symbolic_simplify_false.svg" alt="drawing" width="50%"> 
<img src="Illustrations\figure_chebyshev_MakeFunction_simplify_false.svg" alt="drawing" width="50%"> 
<img src="Illustrations\figure_chebyshev_Exe_simplify_false.svg" alt="drawing" width="50%">



#### Chebyshev benchmarks with simplification on
<img src="Illustrations\figure_chebyshev_Exe_simplify_true.svg" alt="drawing" width="50%">

With simplification on performance of the executable derivative function for Symbolics.jl is slightly better than with simplification off. But simplification processing time is longer.
 
### Spherical Harmonics

The second example is the spherical harmonics function. This is the expression graph for the spherical harmonic function of order 8:
<img src="Illustrations/sphericalharmonics_8.svg" alt="drawing" width="100%">

<details>
    <summary> Source for spherical harmonics benchmark </summary>

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
</details>

As was the case for Chebyshev polynomials the number of paths from the roots to the variables is much greater than the number of nodes in the graph. Once again the y axis is logarithmic.

<img src="Illustrations\figure_spherical_harmonics_Symbolic_simplify_false.svg" alt="drawing" width="50%">
<img src="Illustrations\figure_spherical_harmonics_MakeFunction_simplify_false.svg" alt="drawing" width="50%">
<img src="Illustrations\figure_spherical_harmonics_Exe_simplify_false.svg" alt="drawing" width="50%">
 
 The **Exe** benchmark took many hours to run and was stopped at model size 24 instead of 25 as for the **Symbolic** and **MakeFunction** benchmarks.

</details>


## Limitations
**FD** does not support expressions with conditionals on **FD** variables. For example, you can do this:
```
julia> f(a,b,c) = a< 1.0 ? cos(b) : sin(c)
f (generic function with 2 methods)

julia> f(0.0,x,y)
cos(x)

julia> f(1.0,x,y)
sin(y)
```
but you can't do this:
```
julia> f(a,b) = a < b ? cos(a) : sin(b)
f (generic function with 2 methods)

julia> f(x,y)
ERROR: MethodError: no method matching isless(::FastDifferentiation.Node{Symbol, 0}, ::FastDifferentiation.Node{Symbol, 0})

Closest candidates are:
  isless(::Any, ::DataValues.DataValue{Union{}})
   @ DataValues ~/.julia/packages/DataValues/N7oeL/src/scalar/core.jl:291
  isless(::S, ::DataValues.DataValue{T}) where {S, T}
   @ DataValues ~/.julia/packages/DataValues/N7oeL/src/scalar/core.jl:285
  isless(::DataValues.DataValue{Union{}}, ::Any)
   @ DataValues ~/.julia/packages/DataValues/N7oeL/src/scalar/core.jl:293
  ...
```
This is because the call `f(x,y)` creates an expression graph. At graph creation time the **FD** variables `x,y` are unevaluated variables with no specific value so they cannot be compared with any other value.

The algorithm can be extended to work with conditionals applied to **FD** variables but the processing time and graph size may grow exponentially with conditional nesting depth. A future version may allow for limited conditional nesting. See [Future Work](#FutureWork) for a potential long term solution to this problem.

The preprocessing/compilation step for expressions graphs with more than ≈10⁵ operations may take a minute or more. This is due to two factors. 

The current code is not memory efficient - it allocates much more than necessary which makes it slower than it should be. Future versions will be more memory efficient.

The code uses BitVector for tracking reachability of function roots and variable nodes. This seemed like a good idea when I began and thought **FD** would only be practical for modest size expressions (<10⁴ operations). But, it scaled better than expected and for larger graphs the memory overhead of the BitVector representation is significant. Using Set instead of BitVector for larger graphs should significantly reduce symbolic processing time.


# How it works
The **FD** differentiation algorithm is related to the [D*](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) algorithm but is asymptotically faster so it works on much larger expression graphs. The new algorithms used in **FD** will be described in a soon to be written paper.

**FD** transforms the input expression graph into a derivative graph[^a], and then factors this graph to generate an efficient expression for the derivative. This is fundamentally different from forward and reverse automatic differentiation. 

The efficiency of **FD** comes from analysis of the graph structure of the function rather than sophisticated algebraic simplification rules. By default **FD** applies only these algebraic simplications[^c] to expressions:
* `x×0=>0`
* `x×1=>x`
* `x/1=>x`
* `x+0=>x`
* `c₁×c₂=>c₃` for `c₁,c₂,c₃` constants
* `c₁+c₂=>c₃` for `c₁,c₂,c₃` constants
* `c₁×(c₂×x)` => `(c₁×c₂)×x`  for `c₁,c₂` constants

These rules are generally safe in the sense of obeying IEEE floating point arithmetic rules. However if the runtime value of `x` happens to be NaN or Inf the **FD** expression `x*0` will identically return 0, because it will have been rewritten to 0 by the simplification rules. The expected IEEE result is NaN.


<div id="FutureWork"></div>

## Future work
The **FD** algorithm is fast enough to preprocess expression graphs with ≈10⁵ operations in approximately 1 minute on a modern laptop (processing time scales non-linearly so smaller graphs take much less time).

However, LLVM compile time can be significant at this scale. For expressions this size and larger [DynamicExpressions.jl](https://github.com/SymbolicML/DynamicExpressions.jl) might be a better tradeoff between compile and execution time. This would be especially useful when your function is changing frequently so compilation overhead cannot be amortized across many derivative evaluations.

In its current form **FD** cannot differentiate expressions with conditionals involving **FD** variables. The algorithm can be extended to allow this but symbolic processing can scale exponentially with the nesting depth of conditionals. For small nesting depths this might be acceptable so a future version of FD might support limited nested conditionals. 

However, a better approach might be to use FD as a processing step in a tracing JIT compiler, applying **FD** to the basic blocks detected and compiled by the JIT. These basic blocks do not have branches. Many programs could be differentiated competely automatically by this method. I'm not a compiler expert so it is unlikely I will do this by myself. But contact me if *you* are a compiler expert and want a cool project to work on.

[^c]: More rules may be added in future versions of FD to improve efficiency.

[^b]: I am working with the SciML team to see if it is possible to integrate **FD** differentiation directly into Symbolics.jl.

[^a]: See the [D* ](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) paper for an explanation of derivative graph factorization. 

[^d]: See the Lagrangian dynamics example in the [D* ](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) paper.
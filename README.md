# FastSymbolicDifferentiation

[![Build Status](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package is not yet on the Julia registry and is currently Beta software. I hope to have it registered in a few weeks. Come back then for better documentation, fewer bugs, and more features.

This is a package for computing symbolic derivatives quickly and for generating efficient executables to evaluate those derivatives. It uses a new algorithm, called **FSD**, which is related to the [D* ](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) algorithm but is asymptotically  faster.  

For expression graphs with many common subexpressions (where each node in the expression graph has more than one parent on average) **FSD** may compute symbolic derivatives much more quickly than conventional computer algegra systems such as Symbolics.jl or Mathematica. The generated executables may also be significantly faster. 

If your function is small or tree like (where each node in the expression graph has one parent on average) then Symbolics.jl may outperform **FSD**. There is no simple rule to determine whether **FSD** will do better than Symbolics.jl. Try them both and choose the faster one.

**FSD** can be used standalone if all you need is a derivative or with Symbolics.jl if you need to do further analysis on the symbolic derivative. Converting between Symbolics.jl and **FSD** expression forms is straightforward. However, because of the tree based representation used by Symbolics.jl expression size can grow significantly when converting from **FSD** to Symbolics.jl forms.


Unlike forward and reverse automatic differentiation you don't have to choose which differentiation algorithm to use based on the graph structure. **FSD** automatically generates efficient derivatives for arbitrary function types: ℝ¹->ℝ¹, ℝ¹->ℝᵐ, ℝⁿ->ℝ¹, and ℝⁿ->ℝᵐ, m≠1,n≠1. Its efficiency comes from analysis of the graph structure of the function rather than sophisticated algebraic simplification rules. By default **FSD** applies only these algebraic simplications[^1] to expressions:
* x×0=>0
* x×1=>x
* x/1=>x
* x+0=>x
* c₁×c₂=>c₃ for c₁,c₂,c₃ constants
* c₁+c₂=>c₃ for c₁,c₂,c₃ constants
* c₁×(c₂×x) => (c₁×c₂)×x  for c₁,c₂ constants


These rules are generally safe in the sense of obeying IEEE floating point arithmetic rules. However if the runtime value of x happens to be NaN or Inf the **FSD** expressions x*0 and x+0 will identically return 0, because they will have been rewritten to 0 by the simplification rules. The expected IEEE result in these cases would be NaN. 

## Benchmarks

The FSDBenchmark subdirectory has several benchmark functions you can use to compare the performance of [Symbolics.jl](https://symbolics.juliasymbolics.org/dev/) to FastSymbolicDifferentiation.jl on your computer. There are three types of benchmars:Symbolic,MakeFunction, and Exe. 

The Symbolic benchmark compares the time required to compute just the symbolic form of the derivative. The Symbolic benchmark can be run with simplification turned on or off for Symbolics.jl. If simplification is on then computation time can be extremely long but the resulting expression might be simpler and faster to execute.

The MakeFunction benchmark generates a Julia Expr but does not compile it. This separates the time required to create the function from the time LLVM takes to compile it. As symbolic expressions become large LLVM compile time and memory usage both increase dramatically. For both of the example benchmarks shown below LLVM ran out of memory at relatively small problem sizes.

The Exe function measures just the time required to execute the compiled function using an in-place matrix.

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
The function is memoized for efficiency. Symbolics.jl can simplify the resulting expression graphs to a simple polynomial form when full simplification is turned on. This yields efficient executables but the symbolic processing can take a very long time.

All benchmarks show the ratio of time taken by Symbolics.jl to FastSymbolicDifferentiation.jl. Numbers greater than 1 mean FastSymbolicDifferentiation is faster.

<img src="FSDBenchmark\Data\figure_chebyshev_Symbolic.svg" alt="drawing" width="50%">
<img src="FSDBenchmark\Data\figure_chebyshev_MakeFunction.svg" alt="drawing" width="50%">
<img src="FSDBenchmark\Data\figure_chebyshev_Exe.svg" alt="drawing" width="50%">

<img src="FSDBenchmark\Data\figure_spherical_harmonics_Symbolic.svg" alt="drawing" width="50%">
<img src="FSDBenchmark\Data\figure_spherical_harmonics_MakeFunction.svg" alt="drawing" width="50%">
<img src="FSDBenchmark\Data\figure_spherical_harmonics_Exe.svg" alt="drawing" width="50%">


[^1]: More rules may be added in future versions of FSD to improve efficiency.

# FastSymbolicDifferentiation

[![Build Status](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a package for computing symbolic derivatives quickly and for generating efficient executables to evaluate those derivatives. It uses a new algorithm, called **FSD**, which is related to the [D* ](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) algorithm but is asymptotically  faster[^1].  

For expression graphs with many common subexpressions (where each node in the expression graph has more than one parent on average) **FSD** may compute symbolic derivatives much more quickly than conventional computer algegra systems such as Symbolics.jl or Mathematica. The generated executables may also be significantly faster. 

If your function is small or tree like (where each node in the expression graph has one parent on average) then Symbolics.jl may outperform **FSD**. There is no simple rule to determine whether **FSD** will do better than Symbolics.jl. Try them both and choose the faster one.

**FSD** can be used standalone if all you need is a derivative or with Symbolics.jl if you need to do further analysis on the symbolic derivative. Converting between Symbolics.jl and **FSD** expression forms is straightforward.


Unlike forward and reverse automatic differentiation you don't have to choose which differentiation algorithm to use based on the graph structure. **FSD** automatically generates efficient derivatives for arbitrary function types: ℝ¹->ℝ¹, ℝ¹->ℝᵐ, ℝⁿ->ℝ¹, and ℝⁿ->ℝᵐ, m≠1,n≠1. Its efficiency comes from analysis of the graph structure of the function rather than sophisticated algebraic simplification rules. By default **FSD** applies only these algebraic simplications[^2] to expressions:
* x*0=>0
* x*1=>x
* x/1=>x
* x+0=>x
* c1*c2=>c3 for c₁,c₂,c₃ constants
* c1+c2=>c3 for c₁,c₂,c₃ constants


These rules are generally safe in the sense of obeying IEEE floating point arithmetic rules. However if the runtime value of x happens to be NaN or Inf the **FSD** expressions x*0 and x+0 will identically return 0, because they will have been rewritten to 0 by the simplification rules. The expected IEEE result in these cases would be NaN. 

[^1]: O(m²+n²)|E| for FSD versus O(mn)|E^2| for D* where n is the domain dimension of the ℝⁿ->ℝᵐ function being differentiated, m is the codomain dimension, and |E| is the number of edges in the expression graph. Except for trivial graphs m,n are both typically much smaller than |E|.
[^2]: More rules may be added in future versions of FSD to improve efficiency.

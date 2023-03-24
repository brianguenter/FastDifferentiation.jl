# FastSymbolicDifferentiation

[![Build Status](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a package for computing symbolic derivatives quickly and for generating efficient executables to evaluate those derivatives. It uses a new algorithm, called **FSD**, which is related to the [D* ](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) algorithm but is asymptotically much faster[^1].  

For expression graphs with many common subexpressions **FSD** may compute symbolic derivatives much more quickly than conventional computer algegra systems such as Symbolics.jl or Mathematica. The generated executables may also be significantly faster.

Unlike forward and reverse automatic differentiation the user does not have to choose which differentiation algorithm to use based on the graph structure. **FSD** automatically analyzes the graph structure and efficiently extracts shared derivative terms in the ℝⁿ->ℝᵐ function.

**FSD** is designed to be compatible with Symbolics.jl so it is straightforward to convert Symbolics.jl statements into **FSD** compatible forms and vice versa.

[^1]: O(m+n)|E| for FSD versus O(mn)|E^2| for D* where n is the domain dimension of the ℝⁿ->ℝᵐ function being differentiated, m is the codomain dimension, and |E| is the number of edges in the expression graph.

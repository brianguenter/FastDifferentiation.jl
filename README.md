# FastSymbolicDifferentiation

[![Build Status](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a package for computing symbolic derivatives quickly and for generating efficient executables to evaluate those derivatives. It uses a new algorithm, called FSD, which is related to the [D* ](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) algorithm but is asymptotically much faster: O(m+n)|E| for FSD versus O(mn)|E^2| for D* where n is the domain dimension of the ℝⁿ->ℝᵐ function being differentiated, m is the codomain dimension, and |E| is the number of edges in the expression graph.

It can compute symbolic derivatives much more quickly than conventional computer algegra systems such as Symbolics.jl or Mathematica. The generated executables can also be significantly more efficient.

FSD is neither forward nor reverse automatic differentiation. The FSD algorithm automatically finds shared derivative terms in ℝⁿ->ℝᵐ functions; Unlike forward and reverse automatic differntiation the user does not have to choose which differentiation algorithm to use based on the graph structure. FSD automatically analyzes the graph structure and extracts shared derivative terms in ℝⁿ->ℝᵐ functions

# FastSymbolicDifferentiation

[![Build Status](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a package for computing symbolic derivatives quickly and for generating efficient executables to evaluate those derivatives. It uses a new algorithm, called FSD, which is related to the [D* ](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) algorithm but is asymptotically much faster: O(m+n)|E| for FSD versus O(mn)|E^2| for D* where n is the domain dimension of the ℝⁿ->ℝᵐ being differentiated, m is the codomain dimension, and |E| is the number of edges in the expression graph.

For expression graphs with many common subexpressions it can also be much faster than conventional computer algegra systems such as Symbolics.jl or Mathematica.

FSD is neither forward nor reverse automatic differentiation. The FSD algorithm automatically finds shared derivative terms in ℝⁿ->ℝᵐ functions.

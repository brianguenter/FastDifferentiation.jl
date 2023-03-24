# FastSymbolicDifferentiation

[![Build Status](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brianguenter/FastSymbolicDifferentiation.jl/actions/workflows/CI.yml?query=branch%3Amain)

A new algorithm for computing symbolic derivatives quickly. It is related to the [D* ](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) algorithm but is asymptotically much faster: O(m+n)|E| for FSD versus O(mn)|E^2| for D*. It is also much faster in practice for graphs of all sizes.

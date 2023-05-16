# This script is for setting up the REPL environment for development of the Differentation package. No use in production.

using FastSymbolicDifferentiation
using FastSymbolicDifferentiation.FSDTests
using StaticArrays

using Symbolics: @variables

@variables x y z
nx = Node(x)
ny = Node(y)
nz = Node(z)

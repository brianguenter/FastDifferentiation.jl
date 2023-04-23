using FSDBenchmark
using FastSymbolicDifferentiation
using Symbolics
using StaticArrays
using FSDBenchmark.SimpsonHermite

@variables x y z t
nx = Node(x)
ny = Node(y)
nz = Node(z)
nt = Node(t)
q = function_of(:q, nt)

gr = SiH_test();


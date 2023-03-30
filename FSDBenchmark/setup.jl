using FSDBenchmark
using FastSymbolicDifferentiation
using Symbolics

@variables x y z t
nx = Node(x)
ny = Node(y)
nz = Node(z)
nt = Node(t)
q = function_of(:q, nt)

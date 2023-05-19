using FSDBenchmark
using FastSymbolicDifferentiation
using Symbolics: @variables
using StaticArrays
using FSDBenchmark.SimpsonHermite

@variables x y z t
nx = Node(x)
ny = Node(y)
nz = Node(z)
nt = Node(t)

# gr = SiH_test();
# time_test(deepcopy(gr))


#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using FastSymbolicDifferentiation.TestCases
using StaticArrays
using FiniteDifferences
using Profile
using PProf
using .TestCases

function test()
    order = 7
    Symbolics.@variables x y z

    derivs = SHDerivatives(order, x, y, z)
end
export test

function profile()
    graph, qx, qy, qz = to_graph(25)
    Profile.clear()
    @profile symbolic_jacobian!(graph, [Node(qx), Node(qy), Node(qz)])
    pprof()
end
export profile

#changed
#change
export test

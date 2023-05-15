#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    @variables x
    nx = Node(x)
    zr = Node(1.0)

    graph = DerivativeGraph([nx, zr])
    jac = symbolic_jacobian!(graph, [nx])

end



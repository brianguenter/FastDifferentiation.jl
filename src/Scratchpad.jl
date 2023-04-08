#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables x y

    nx = Node(x)
    ny = Node(y)
    n2 = nx * ny
    n4 = n2 * ny
    n5 = n2 * n4

    graph = DerivativeGraph([n4, n5])

    df21(x, y) = 2 * x * y^3
    df22(x, y) = 3 * x^2 * y^2
    df11(x, y) = y^2
    df12(x, y) = 2 * x * y
    factor!(graph)
    display(symbolic_jacobian!(graph, [nx, ny]))
    Vis.draw_dot(graph, reachability_labels=false)

end
export test

#changed
#change
export test

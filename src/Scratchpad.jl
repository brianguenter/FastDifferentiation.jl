#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
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
    df22(x, y) = 4 * x^2 * y^2
    df11(x, y) = y^2
    df12(x, y) = 2 * x * y

    correct_jacobian = [df11 df12; df21 df22]
    copy_jac = _symbolic_jacobian(graph, [nx, ny])
    jac = _symbolic_jacobian!(graph, [nx, ny])

    @assert all(copy_jac .== jac) #make sure the jacobian computed by copying the graph has the same variables as the one computed by destructively modifying the graph

    computed_jacobian = make_function(jac, [nx, ny])

    #verify the computed and hand caluclated jacobians agree.
    for _x in -1.0:0.01:1.0
        for _y in -1.0:0.3:1.0
            for index in CartesianIndices(correct_jacobian)
                @assert isapprox(correct_jacobian[index](_x, _y), computed_jacobian(_x, _y)[index])
            end
        end
    end
end




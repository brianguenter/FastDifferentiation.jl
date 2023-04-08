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

    correct_jacobian = [df11 df12; df21 df22]

    symbolic = symbolic_jacobian!(graph, [nx, ny])
    computed_jacobian = make_function.(symbolic, Ref([x, y]))

    println(symbolic[2, 2])
    #verify the computed and hand caluclated jacobians agree.
    for x in -1.0:0.5:1.0
        for y in -1.0:0.5:1.0
            for index in CartesianIndices(correct_jacobian)
                # println("correct $(correct_jacobian[index](x,y)) computed $(computed_jacobian[index](x,y))")
                # @assert isapprox(correct_jacobian[index](x, y), computed_jacobian[index](x, y))
            end
        end
    end

end
export test

#changed
#change
export test

#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    fsd_graph = FSD_spherical_harmonics(10)
    fsd_func = make_function(fsd_graph, variables(fsd_graph))

    #hand computed derivative for order = 3
    # correct_derivatives(x, y, z) = [
    #     0.0 0.0 0.0
    #     0.0 0.0 0.0
    #     0.0 0.0 1.422074017360395
    #     0.0 0.0 0.0
    #     0.0 0.0 0.0
    #     (-0.3642161365141257*3*z) 0.0 (-0.3642161365141257*3*x)
    #     0.0 0.0 (2.0197963935867267*3*z)
    #     0.0 (-0.3642161365141257*3*z) (-0.3642161365141257*3*y)
    #     0.0 0.0 0.0
    # ]

    sym_func = jacobian_function!(fsd_graph, variables(fsd_graph))
end
export test



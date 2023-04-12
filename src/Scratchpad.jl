#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()

    fsd_graph, x, y, z = to_graph(10)
    # fsd_func = make_function(fsd_graph, Node.([x, y, z]))

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

    sym = symbolic_jacobian!(fsd_graph)
end
export test


function make_subdata()
    @variables x y

    nx = Node(x)
    ny = Node(y)
    n2 = nx * ny
    n4 = n2 * ny
    n5 = n2 * n4


    graph = DerivativeGraph([n4, n5])
    subs, _ = compute_factorable_subgraphs(graph)

    _5_3 = subs[1]
    @assert (5, 3) == vertices(_5_3)
    e5_3 = edges(graph, 5, 3)[1]
    e3_4 = edges(graph, 3, 4)[1]

    return path(_5_3, e3_4)
end
export make_subdata

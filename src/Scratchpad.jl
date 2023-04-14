#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()

    fsd_graph, x, y, z = to_graph(5)

    factor!(fsd_graph)
    Vis.write_dot("factored.pdf", fsd_graph)
    # sym = symbolic_jacobian!(fsd_graph)
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

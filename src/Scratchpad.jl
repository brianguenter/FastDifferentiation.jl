#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    @variables x
    nx = Node(x)
    println(Chebyshev(2, nx))
    fsd_graph = chebyshev_graph(2)
    Vis.draw_dot(fsd_graph)


    sym_func = symbolic_jacobian!(fsd_graph)
end
export test

function test2()
    @variables x

    nx = Node(x)
    n2 = cos(nx)
    n3 = n2 * n2

    gr = DerivativeGraph(n3)
    Vis.draw_dot(gr)
    sub = postdominator_subgraph(gr, 1, 3, BitVector([1]), BitVector([1]), BitVector([1]))
    evaluate_branching_subgraph(sub)
end
export test2


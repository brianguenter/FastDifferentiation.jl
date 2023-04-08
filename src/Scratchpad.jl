#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables x y

    nx1 = Node(x)
    ny2 = Node(y)
    n3 = nx1 * ny2
    n4 = n3 * ny2
    n5 = n3 * n4


    graph = DerivativeGraph([n5, n4])
    tmp = postdominator_subgraph(graph, 2, 4, BitVector([0, 1]), BitVector([0, 1]), BitVector([0, 1]))
    println(summarize(tmp))
    factor_subgraph!(tmp)
    Vis.draw_dot(graph)
    # @test length(edges(graph, 2, 4)) == 1
end
export test

#changed
#change
export test

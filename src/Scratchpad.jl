#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables x y

    q = function_of(:q, x, y)
    f = x * q + y * q
    graph = DerivativeGraph([f])
    factor!(graph)
    remove_dangling_edges!(graph)
    Vis.draw(graph)
    value()
    # symbolic_jacobian!(graph)
end
export test

#changed
#change
export test

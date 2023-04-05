#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    _, graph, _, _ = simple_dominator_graph()
    factor!(graph)
    fedge = edges(graph, 1, 4)[1]
    dfsimp = make_function(value(fedge))
    _, graph, _, _ = simple_dominator_graph()
    origfsimp = make_function(root(graph, 1))
    @assert isapprox(central_fdm(5, 1)(origfsimp, 3), dfsimp(3))

    graph = complex_dominator_graph()
    factor!(graph)
    fedge = edges(graph, 1, 8)[1]
    df = make_function(value(fedge))

    graph = complex_dominator_graph()
    origf = make_function(root(graph, 1))
    @assert isapprox(central_fdm(5, 1)(origf, 3), df(3))
end
export test

#changed
#change
export test

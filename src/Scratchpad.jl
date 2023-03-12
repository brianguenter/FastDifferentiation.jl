#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using .TestCases
using StaticArrays
using .SphericalHarmonics
using FiniteDifferences
using Profile
using PProf

function test()
    _, graph, _, _ = simple_dominator_graph()
    factor!(graph)
    fedge = edges(graph, 1, 4)[1]
    dfsimp = dag_to_function(edge_value(fedge))
    _, graph, _, _ = simple_dominator_graph()
    origfsimp = dag_to_function(root(graph, 1))
    @assert isapprox(central_fdm(5, 1)(origfsimp, 3), dfsimp(3))

    graph = complex_dominator_graph()
    factor!(graph)
    fedge = edges(graph, 1, 8)[1]
    df = dag_to_function(edge_value(fedge))

    graph = complex_dominator_graph()
    origf = dag_to_function(root(graph, 1))
    @assert isapprox(central_fdm(5, 1)(origf, 3), df(3))
end
export test

function profile()
    graph, qx, qy, qz = to_graph(25)
    Profile.clear()
    @profile symbolic_jacobian!(graph, [Node(qx), Node(qy), Node(qz)])
    pprof()
end
export profile

export test

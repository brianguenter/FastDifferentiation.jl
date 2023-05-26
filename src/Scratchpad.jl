#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    x, graph, _, _ = simple_dominator_graph()
    nx = Node(x)
    factor!(graph)
    fedge = edges(graph, 1, 4)[1]
    tmp0 = make_function([value(fedge)], [nx])
    dfsimp(x) = tmp0(x)[1]
    x, graph, _, _ = simple_dominator_graph() #x is a new variable so have to make a new Node(x)
    nx = Node(x)
    tmp00 = make_function([root(graph, 1)], [nx])
    origfsimp(x) = tmp00(x)[1]
    @assert isapprox(central_fdm(5, 1)(origfsimp, 3), dfsimp(3)[1])

    graph = complex_dominator_graph()
    factor!(graph)
    fedge = edges(graph, 1, 8)[1]
    tmp1 = make_function([value(fedge)], variables(graph))
    df(x) = tmp1(x)[1]

    graph = complex_dominator_graph()
    tmp2 = make_function([root(graph, 1)], variables(graph))
    origf(x) = tmp2(x)[1]

    for test_val in -3.0:0.013:3.0
        @assert isapprox(central_fdm(5, 1)(origf, test_val), df(test_val)[1])
    end
end




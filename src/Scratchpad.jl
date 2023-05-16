#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    @variables x, y

    #f = [cos(x) * sin(y), cos(x) + sin(y)]
    nx = Node(x)
    ny = Node(y)
    cosx = Node(cos, nx)
    sinx = Node(sin, nx)
    partial_cosx = Node(-, sinx)
    siny = Node(sin, ny)
    partial_siny = Node(cos, ny)
    ctimess = cosx * siny
    partial_times = [siny, cosx]
    cpluss = cosx + siny
    partial_plus = [Node(1), Node(1)]
    roots = [ctimess, cpluss]
    grnodes = [nx, ny, cosx, siny, cpluss, ctimess]

    correct_postorder_numbers = Dict((nx => 1, cosx => 2, ny => 3, siny => 4, ctimess => 5, cpluss => 6))

    graph = DerivativeGraph(roots)

    @assert all([correct_postorder_numbers[node] == postorder_number(graph, node) for node in grnodes])

    correct_partials = Dict((cosx => [partial_cosx], siny => [partial_siny], ctimess => partial_times, cpluss => partial_plus))
end
export test



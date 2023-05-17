#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    @variables x, y, z
    fsd_graph = spherical_harmonics(FastSymbolic(), 25, x, y, z)
    println("roots $(length(roots(fsd_graph)))")
    println(number_of_operations(fsd_graph))

    sparse = sparse_symbolic_jacobian!(fsd_graph, variables(fsd_graph))
    println(number_of_operations(sparse))
end
export test



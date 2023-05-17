#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    @variables x, y, z
    fsd_graph = spherical_harmonics(FastSymbolic(), 10, x, y, z)
    sprse = sparse_symbolic_jacobian!(fsd_graph, variables(fsd_graph))
end
export test



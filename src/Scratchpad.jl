#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    fsd_graph = FSD_spherical_harmonics(10)
    return sparse_symbolic_jacobian!(fsd_graph, variables(fsd_graph))
end
export test



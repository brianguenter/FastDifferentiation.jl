#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    fsd_graph = spherical_harmonics(FastSymbolic(), 10)
    fsd_func = make_function(fsd_graph, variables(fsd_graph))

    sym_func = _jacobian_function!(fsd_graph, variables(fsd_graph))
end




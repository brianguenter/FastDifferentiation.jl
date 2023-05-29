#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    fsd_graph = spherical_harmonics(FastSymbolic(), 10)
    mn_func = make_function(roots(fsd_graph), variables(fsd_graph))
    fsd_func(vars...) = vec(mn_func(vars))

    graph_vars = variables(fsd_graph)

    sym_func = make_function(symbolic_jacobian(roots(fsd_graph), graph_vars),graph_vars)
end




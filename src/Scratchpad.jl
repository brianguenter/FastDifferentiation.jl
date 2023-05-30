#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences



function test()
    order = 10

    fsd_graph = spherical_harmonics(FastSymbolic(), order)
    fsd_func = roots(fsd_graph)
    func_vars = variables(fsd_graph)

    Jv, v_vars = jacobian_times_v(fsd_func, func_vars)

end




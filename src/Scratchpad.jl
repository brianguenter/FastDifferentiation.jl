#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    chebyshev_order = 20
    fsd_graph = FSD_chebyshev(chebyshev_order)
    fsd_func = make_function(fsd_graph)

    func_wrap(x) = fsd_func(x)[1]

    sym_func = jacobian_function!(fsd_graph, in_place=false)

    for xr in -1.0:0.214:1.0
        finite_diff = central_fdm(12, 1, adapt=3)(func_wrap, xr)

        symbolic = sym_func(xr)

        @assert isapprox(symbolic[1, 1], finite_diff[1], rtol=1e-8)
    end
end
export test



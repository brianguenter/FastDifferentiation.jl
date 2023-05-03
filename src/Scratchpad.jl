#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    fsd_graph = chebyshev_graph(5)
    fsd_func = make_function(fsd_graph)

    func_wrap(x) = fsd_func(x)[1]

    sym_func = jacobian_function!(fsd_graph)

    for xr in -1.0:0.214:1.0
        finite_diff = central_fdm(12, 1, adapt=3)(func_wrap, xr)

        symbolic = sym_func(xr)

        @assert isapprox(symbolic[1, 1], finite_diff[1], rtol=1e-8)
    end

    tmp = Matrix{Float64}(undef, 1, 1)
    fsd_graph = chebyshev_graph(20)
    sym_func = jacobian_function!(fsd_graph)
end
export test

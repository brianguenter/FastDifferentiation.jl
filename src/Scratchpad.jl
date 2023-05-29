#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    order = 10

    fsd_graph = spherical_harmonics(FastSymbolic(), order)
    fsd_func = roots(fsd_graph)
    func_vars = variables(fsd_graph)

    Jᵀv, r_vars = jacobian_transpose_v(fsd_func, func_vars)

    Jᵀv_slow = convert.(Node, transpose(symbolic_jacobian(fsd_func, func_vars)) * r_vars)
    both_vars = [func_vars; r_vars]
    slow_symbolic = reshape(Jᵀv_slow, (length(Jᵀv_slow), 1))

    slow = make_function(slow_symbolic, both_vars)
    fast = make_function(reshape(Jᵀv, (length(Jᵀv), 1)), both_vars)

    for _ in 1:100
        input = rand(length(func_vars) + length(r_vars))
        slow_val = slow(input)
        fast_val = fast(input)

        @assert isapprox(slow_val, fast_val, rtol=1e-8)
    end

    fast2 = jacobian_transpose_v_exe(fsd_func, func_vars)

    for _ in 1:100
        xin = rand(length(fsd_func))
        vin = rand(codomain_dimension(fsd_graph))
        slow_val = slow([xin; vin])
        fast_val = fast2([xin; vin])

        @assert isapprox(slow_val, fast_val, rtol=1e-8)
    end
end




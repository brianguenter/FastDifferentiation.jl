#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    order = 30

    fsd_graph = spherical_harmonics(FastSymbolic(), order)
    fsd_func = roots(fsd_graph)
    func_vars = variables(fsd_graph)

    Jv, v_vars = jacobian_times_v(fsd_func, func_vars)

    #compute the product the slow way
    Jv_slow = convert.(Node, symbolic_jacobian(fsd_func, func_vars) * v_vars)
    both_vars = [func_vars; v_vars]
    slow_symbolic = reshape(Jv_slow, (length(Jv_slow), 1))

    slow = make_function(slow_symbolic, both_vars)
    fast = make_function(reshape(Jv, (length(Jv), 1)), both_vars)

    for _ in 1:100
        input = rand(length(func_vars) + length(v_vars))
        slow_val = slow(input...)
        fast_val = fast(input...)

        @assert isapprox(slow_val, fast_val, rtol=1e-8) "slow_val $slow_val \n fast_val $fast_val"
    end

    println("func ops = $(number_of_operations(fsd_graph))")
    println("jac ops = $(number_of_operations(symbolic_jacobian(roots(fsd_graph),variables(fsd_graph))))")
    println("jtv ops = $(number_of_operations(Jv))")

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
        slow_val = slow(input...)
        fast_val = fast(input...)

        @assert isapprox(slow_val, fast_val, rtol=1e-8) "slow_val $slow_val \n fast_val $fast_val"
    end
end




#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    order = 10

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

        @assert isapprox(slow_val, fast_val, rtol=1e-9) "slow_val $slow_val \n fast_val $fast_val"
    end

    fast2 = jacobian_times_v_exe(fsd_func, func_vars)

    for _ in 1:100
        xin = rand(length(fsd_func))
        vin = rand(domain_dimension(fsd_graph))
        slow_val = slow([xin; vin]...)
        fast_val = fast2(xin, vin)

        @assert isapprox(slow_val, fast_val, rtol=1e-8) "slow_val $slow_val \n fast_val $fast_val"
    end
end




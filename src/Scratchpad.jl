#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    # order = 4

    # fsd_graph = spherical_harmonics(FastSymbolic(), order)
    # fsd_func = roots(fsd_graph)
    # func_vars = variables(fsd_graph)

    @variables x y
    nx, ny = Node.((x, y))
    fsd_func = [nx^2 * ny, nx * ny^2]
    func_vars = [nx, ny]

    Jᵀv, v_vars = jacobian_transpose_v(fsd_func, func_vars)

    #compute the product the slow way
    Jᵀv_slow = convert.(Node, symbolic_jacobian(fsd_func, func_vars) * v_vars)
    both_vars = [func_vars; v_vars]
    slow_symbolic = reshape(Jᵀv_slow, (length(Jᵀv_slow), 1))
    println("correct\n")
    display(slow_symbolic)
    println("\nunder test\n")
    display(Jᵀv)
    slow = make_function(slow_symbolic, both_vars)
    fast = make_function(reshape(Jᵀv, (length(Jᵀv), 1)), both_vars)

    println(methods(slow))
    for _ in 1:100
        input = rand(length(func_vars) + length(v_vars))
        slow_val = slow(input...)
        fast_val = fast(input...)

        @assert isapprox(slow_val, fast_val, rtol=1e-10) "slow_val $slow_val \n fast_val $fast_val"
    end
end




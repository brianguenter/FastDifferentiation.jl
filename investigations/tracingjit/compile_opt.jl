# Does a lower LLVM optimization level cut the compile time of the big generated Jacobian functions, and what does it cost
# at run time?  The same generated expression is compiled (a) the default way with RuntimeGeneratedFunctions and (b) by
# evaluating it in a module whose optimization level is set with Base.Experimental.@optlevel.
#   julia --project=<env> investigations/tracingjit/compile_opt.jl 4

include(joinpath(@__DIR__, "prototype.jl"))
using .TracingPrototype
using Random, Printf, LinearAlgebra

const TP = TracingPrototype
rng = MersenneTwister(1)
n = parse(Int, isempty(ARGS) ? "4" : ARGS[1])

function jacobi_svdvals(A; tol=1e-13, maxsweeps=60)
    m, n = size(A)
    B = copy(A)
    for sweep in 1:maxsweeps
        rotated = false
        for p in 1:n-1, q in p+1:n
            α = zero(eltype(B)); β = zero(eltype(B)); γ = zero(eltype(B))
            for i in 1:m
                α += B[i, p] * B[i, p]; β += B[i, q] * B[i, q]; γ += B[i, p] * B[i, q]
            end
            if abs(γ) > tol * sqrt(α * β)
                θ = atan(2γ, β - α) / 2
                c = cos(θ); s = sin(θ)
                for i in 1:m
                    bp = B[i, p]; bq = B[i, q]
                    B[i, p] = c * bp - s * bq
                    B[i, q] = s * bp + c * bq
                end
                rotated = true
            end
        end
        rotated || break
    end
    return sort([sqrt(sum(abs2, B[:, j])) for j in 1:n]; rev=true)
end

function timeit(thunk)
    thunk(); thunk()
    k = 1
    while true
        t = @elapsed for _ in 1:k
            thunk()
        end
        t > 0.2 && return t / k
        k *= 4
    end
end
runner(f, res, x, reps) = (for _ in 1:reps
    f(res, x)
end; nothing)

A = randn(rng, n, n) + n * I
f = A -> sum(jacobi_svdvals(A))
tj = TracingJacobian(f; in_place=true, guard_mode=:replay) # no guards: the generated function is only the Jacobian
res = zeros(1, n * n)
tj(res, A)
tr = tj.traces[1]
x = TP.flatten((A,))
println("nuclear norm gradient, $(n)x$(n): ", TP.trace_report(tr))
@printf("RuntimeGeneratedFunction (default)        compile + first call %8.2f s   call %9.2f us\n", tr.stats.first_call_s, 1e6 * timeit(() -> tr.fn(res, x)))

expr, ns, ng = TP.trace_function_expr(Pair{TP.Node,Bool}[], tr.J, tj.vars; in_place=true)
for L in (2, 1, 0)
    m = Module(Symbol("Gen", L))
    Core.eval(m, :(Base.Experimental.@optlevel $L))
    res2 = zeros(1, n * n)
    t0 = time_ns()
    fn = Core.eval(m, expr)
    Base.invokelatest(fn, res2, x) # first call: compiles
    tc = (time_ns() - t0) / 1e9
    reps = 2000
    Base.invokelatest(runner, fn, res2, x, 3)
    t = @elapsed Base.invokelatest(runner, fn, res2, x, reps)
    @printf("eval in a module with @optlevel %d          compile + first call %8.2f s   call %9.2f us   (max difference from the default %.1e)\n", L, tc, 1e6 * t / reps, maximum(abs.(res2 .- res)))
end

# ---- the same statements compiled in chunks, values that cross chunks go through a workspace vector
for cs in (500, 1000, 2000, 5000)
    t0 = time_ns()
    cf, nst, ng = TP.chunked_trace_function(Pair{TP.Node,Bool}[], tr.J, tj.vars; chunk_size=cs)
    t_build = (time_ns() - t0) / 1e9
    res3 = zeros(1, n * n)
    t0 = time_ns()
    cf(res3, x) # first call: compiles every chunk
    tc = (time_ns() - t0) / 1e9
    t = timeit(() -> cf(res3, x))
    @printf("chunked, %5d statements per chunk (%3d chunks, workspace %6d)  compile + first call %8.2f s   call %9.2f us   (max difference %.1e)\n",
        cs, length(cf.chunks), length(cf.w), tc, 1e6 * t, maximum(abs.(res3 .- res)))
end

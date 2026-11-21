# Compile time of a big trace when it is split into chunks. Does not compile the monolithic version.
#   julia --project=<env> investigations/tracingjit/compile_chunks.jl 6 5000 10000
include(joinpath(@__DIR__, "prototype.jl"))
using .TracingPrototype
using Random, Printf, LinearAlgebra
const TP = TracingPrototype
rng = MersenneTwister(1)
n = parse(Int, ARGS[1])
sizes = parse.(Int, ARGS[2:end])

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
timeit(thunk) = (thunk(); thunk(); k = 1; while true
    t = @elapsed for _ in 1:k
        thunk()
    end
    t > 0.2 && return t / k
    k *= 4
end)

A = randn(rng, n, n) + n * I
f = A -> sum(jacobi_svdvals(A))
tj = TracingJacobian(f; in_place=true, guard_mode=:replay)
x = TP.flatten((A,))
TP.init!(tj, (A,))
tr = TP.record!(tj, TP.unflatten(tj, x), x) # records and differentiates; compiles nothing
println("nuclear norm gradient, $(n)x$(n): ", TP.trace_report(tr))
first_result = nothing
for cs in sizes
    cf, nst, ng = TP.chunked_trace_function(Pair{TP.Node,Bool}[], tr.J, tj.vars; chunk_size=cs)
    res = zeros(1, n * n)
    t0 = time_ns()
    cf(res, x)
    tc = (time_ns() - t0) / 1e9
    t = timeit(() -> cf(res, x))
    global first_result
    first_result === nothing && (first_result = copy(res))
    @printf("chunks of %6d statements (%3d chunks)  compile + first call %8.2f s   call %9.2f us   (difference from the first: %.1e)\n", cs, length(cf.chunks), tc, 1e6 * t, maximum(abs.(res .- first_result)))
end

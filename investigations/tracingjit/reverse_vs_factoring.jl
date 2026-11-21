# Is FD's factored gradient really much bigger than reverse mode on the recorded trace?
# Records the trace of  A -> sum(jacobi_svdvals(A))  and compares, on the identical expression graph:
#   - the trace itself
#   - FD.jacobian (factoring)
#   - FD.reverse_AD (plain reverse accumulation, src/Reverse.jl)
# "ops" = nodes with at least one child, i.e. the statements of generated code (shared subexpressions counted once).
#   julia --project=<env> investigations/tracingjit/reverse_vs_factoring.jl 3 4 5 6

include(joinpath(@__DIR__, "prototype.jl"))
using .TracingPrototype
using FastDifferentiation
using Random, Printf, LinearAlgebra

const TP = TracingPrototype
const FD = FastDifferentiation
const Node = FD.Node

function jacobi_svdvals(A; tol=1e-13, maxsweeps=60)
    m, n = size(A)
    B = copy(A)
    for sweep in 1:maxsweeps
        rotated = false
        for p in 1:n-1, q in p+1:n
            α = zero(eltype(B)); β = zero(eltype(B)); γ = zero(eltype(B))
            for i in 1:m
                α += B[i, p] * B[i, p]
                β += B[i, q] * B[i, q]
                γ += B[i, p] * B[i, q]
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

function record(f, A)
    vars = [Node(Symbol("x", i)) for i in 1:length(A)]
    tracked = reshape([TP.Tracked(A[i], vars[i]) for i in 1:length(A)], size(A))
    rec = TP.Recording()
    TP.CURRENT[] = rec
    out = try
        f(tracked)
    finally
        TP.CURRENT[] = nothing
    end
    return TP.as_nodes(out), vars, rec
end

ops(x) = count(FD.is_tree, FD.all_nodes(x))
allnodes(x) = length(FD.all_nodes(x))

rng = MersenneTwister(7)
sizes = isempty(ARGS) ? [3, 4, 5, 6] : parse.(Int, ARGS)
@printf("%-5s %10s %12s %14s %14s %9s %9s\n", "n", "trace ops", "guards", "FD.jacobian ops", "reverse_AD ops", "jac/trace", "rev/trace")
for n in sizes
    A = randn(rng, n, n) + n * I
    outputs, vars, rec = record(A -> sum(jacobi_svdvals(A)), A)
    t_ops = ops(outputs)
    t_fd = @elapsed J = FD.jacobian(outputs, vars)
    fd_ops = ops(J)
    t_rev = @elapsed R = FD.reverse_AD(outputs[1], vars)
    rev_ops = ops(R)
    @printf("%-5d %10d %12d %14d %14d %9.1f %9.1f    (FD.jacobian %.1f s, reverse_AD %.1f s)\n", n, t_ops, length(rec.guards), fd_ops, rev_ops, fd_ops / t_ops, rev_ops / t_ops, t_fd, t_rev)
end

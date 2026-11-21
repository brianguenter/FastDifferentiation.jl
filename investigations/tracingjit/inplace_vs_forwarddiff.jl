# Fair comparison: the compiled trace evaluated IN PLACE into a preallocated Jacobian (no allocation on a cache hit)
# against ForwardDiff.jacobian! with a prebuilt JacobianConfig and a preallocated Jacobian.
#   julia --project=<env with ForwardDiff, QuadGK, StaticArrays> investigations/tracingjit/inplace_vs_forwarddiff.jl
# The traced function itself may allocate internally when ForwardDiff calls it (matrices, vectors); that is the
# function's own cost, which ForwardDiff cannot avoid and the compiled trace does not pay at all.

include(joinpath(@__DIR__, "prototype.jl"))
using .TracingPrototype
using ForwardDiff, QuadGK, StaticArrays
using Random, Printf, LinearAlgebra

const TP = TracingPrototype
rng = MersenneTwister(11)

function timeit(thunk)
    thunk(); thunk()
    n = 1
    while true
        t = @elapsed for _ in 1:n
            thunk()
        end
        t > 0.3 && return t / n
        n *= 4
    end
end

# f as a function of one flat vector, for ForwardDiff
function flat_function(f, args)
    sizes = [a isa Number ? nothing : size(a) for a in args]
    ranges = TP.arg_ranges(args)
    return function (z)
        newargs = ntuple(i -> sizes[i] === nothing ? z[first(ranges[i])] : reshape(z[ranges[i]], sizes[i]), length(args))
        out = f(newargs...)
        return out isa Number ? [out] : collect(vec(out))
    end
end

function bench(name, f, args)
    tj_ip = TracingJacobian(f; in_place=true)
    tj_op = TracingJacobian(f)
    J_op = tj_op(args...)
    result = similar(J_op)
    tj_ip(result, args...)
    x = TP.flatten(args)
    maxdiff = maximum(abs.(result .- J_op))

    g = flat_function(f, args)
    y = g(x)
    g!(out, z) = (out .= g(z); nothing)
    cfg = ForwardDiff.JacobianConfig(g!, y, x)
    Jf = similar(result)
    ForwardDiff.jacobian!(Jf, g!, y, x, cfg)
    fd_err = maximum(abs.(Jf .- result) ./ max.(abs.(Jf), 1.0))

    eval_jacobian!(tj_ip, result, x)
    alloc_ip = @allocated eval_jacobian!(tj_ip, result, x)
    alloc_op = @allocated eval_jacobian(tj_op, x)
    alloc_fd = @allocated ForwardDiff.jacobian!(Jf, g!, y, x, cfg)

    t_ip = timeit(() -> eval_jacobian!(tj_ip, result, x))
    t_op = timeit(() -> eval_jacobian(tj_op, x))
    t_fd = timeit(() -> ForwardDiff.jacobian!(Jf, g!, y, x, cfg))
    guards = tj_ip.traces[1].stats.guards
    @printf("%-34s %dx%-4d guards %4d | in place %8.2f us %6d B | allocating %8.2f us %6d B | ForwardDiff.jacobian! %8.2f us %6d B | in place is %5.1fx vs ForwardDiff   (in place vs allocating: max diff %.0e, vs ForwardDiff %.1e)\n",
        name, size(result)..., guards, 1e6 * t_ip, alloc_ip, 1e6 * t_op, alloc_op, 1e6 * t_fd, alloc_fd, t_fd / t_ip, maxdiff, fd_err)
end

randm(n) = randn(rng, n, n) + n * I

function f4(u, p)
    for i in 1:30
        u = u > 1.5 ? log(u) : u^2 + p
    end
    return u
end
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

bench("f4 (30 conditionals, 2 in, 1 out)", f4, (7.3, 0.4))
bench("A \\ b, LU 6x6", (A, b) -> A \ b, (randm(6), randn(rng, 6)))
bench("det and inv, 6x6", A -> vcat(det(A), vec(inv(A))), (randm(6),))
bench("cholesky factor, 8x8", A -> vec(Matrix(cholesky(Symmetric(A' * A + I)).L)), (randn(rng, 8, 8),))
bench("Jacobi SVD 4x4 (singular values)", jacobi_svdvals, (randm(4),))
bench("QuadGK cos(a x) on [0,5]", a -> quadgk(x -> cos(a * x), 0.0, 5.0)[1], (20.0,))
bench("StaticArrays eigvals of sym 3x3", a -> collect(eigvals(Symmetric(SMatrix{3,3}(a...)))), (vec(Matrix(Symmetric(randn(rng, 3, 3)))),))

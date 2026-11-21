# Two ways to check the guards of a cached trace:
#   :compiled  the guard conditions are compiled into the trace (the trace recomputes them from the inputs)
#   :replay    the recorded decisions are kept as a packed BitVector (plus the integers produced by Int/floor/...), and f is
#              run again on Probe numbers that compare every decision with the recorded one and stop at the first difference.
#              The trace itself then contains only the Jacobian.
#   julia --project=<env with ForwardDiff, QuadGK, StaticArrays> investigations/tracingjit/replay_vs_compiled.jl

include(joinpath(@__DIR__, "prototype.jl"))
using .TracingPrototype
using ForwardDiff, QuadGK, StaticArrays
using Random, Printf, LinearAlgebra

const TP = TracingPrototype
rng = MersenneTwister(5)

# ============================================================ 1. both modes must give the same Jacobians on every branch
println("1. correctness: the same Jacobians from both modes at random points, which hit many different branches")
bank = [
    ("a < 0 ? x^2 : x^3", (a, x) -> a < 0 ? x^2 : x^3),
    ("min(3t-2, t^2)", (a, t) -> min(3t - 2, t^2)),
    ("max(t, a)^2", (a, t) -> max(t, a)^2),
    ("clamp(t, -a, a) * t", (a, t) -> clamp(t, -a, a) * t),
    ("nested && / ||", (a, t) -> (t < 0 && a > 1) ? t^2 : (t > 1 || a < 0) ? exp(t) : t * a),
    ("30 conditionals in a loop", (u, p) -> (for i in 1:30
        u = u > 1.5 ? log(u) : u^2 + p
    end; u)),
    ("Newton while loop", (a, s0) -> (s = s0 + 0 * a; while abs(s * s - a) > 1e-12
        s = (s + a / s) / 2
    end; s)),
    ("floor(Int, u): B-spline index", (u, p) -> (i = floor(Int, u); (1 - (u - i))^3 / 6 * p[i+1] + (u - i)^3 / 6 * p[i+2])),
    ("Int(mod(u,4)) as an index", (u, x) -> sum(x[i] for i in Int(mod(u, 4))*2+1:Int(mod(u, 4))*2+2)),
]
args_for(name) = occursin("B-spline", name) ? () -> (rand(rng) * 8, randn(rng, 10)) :
                 occursin("mod(u,4)", name) ? () -> (float(rand(rng, 1:12)), randn(rng, 8)) :
                 occursin("30", name) ? () -> (rand(rng) * 40 + 0.1, rand(rng)) :
                 occursin("Newton", name) ? () -> (exp(rand(rng) * 8), 1.0) :
                 () -> (rand(rng) * 4 - 2, rand(rng) * 6 - 3)
for (name, f) in bank
    tc = TracingJacobian(f; max_traces=500, guard_mode=:compiled)
    tr = TracingJacobian(f; max_traces=500, guard_mode=:replay)
    ta = TracingJacobian(f; max_traces=500) # :auto, which may switch from compiled guards to replay part way through
    mk = args_for(name)
    worst = 0.0
    worst_auto = 0.0
    for _ in 1:60
        args = mk()
        Jc = tc(args...)
        Jr = tr(args...)
        Ja = ta(args...)
        worst = max(worst, maximum(abs.(Jc .- Jr)))
        worst_auto = max(worst_auto, maximum(abs.(Jc .- Ja)))
    end
    @printf("   %-30s max |compiled - replay| = %.1e, max |compiled - auto| = %.1e\n   %-30s compiled: %s\n   %-30s replay:   %s\n   %-30s auto:     %s\n",
        name, worst, worst_auto, "", summary_line(tc), "", summary_line(tr), "", summary_line(ta))
end

# ============================================================ 2. speed, in place, against ForwardDiff.jacobian!
println("\n2. speed per call on a cache hit (in place, no allocation of the Jacobian)")
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
    mk(mode) = (tj = TracingJacobian(f; in_place=true, guard_mode=mode); J = tj(similar(zeros(1, 1)), args...) ; tj)
    x = TP.flatten(args)
    tj_c = TracingJacobian(f; in_place=true, guard_mode=:compiled)
    tj_r = TracingJacobian(f; in_place=true, guard_mode=:replay)
    tj_a = TracingJacobian(f; in_place=true) # :auto
    J0 = TracingJacobian(f; guard_mode=:compiled)(args...)
    res_c = similar(J0); res_r = similar(J0); res_a = similar(J0)
    tj_c(res_c, args...)
    tj_r(res_r, args...)
    tj_a(res_a, args...)
    diff = maximum(abs.(res_c .- res_r))
    g = flat_function(f, args)
    y = g(x)
    g!(out, z) = (out .= g(z); nothing)
    cfg = ForwardDiff.JacobianConfig(g!, y, x)
    Jf = similar(res_c)
    ForwardDiff.jacobian!(Jf, g!, y, x, cfg)
    eval_jacobian!(tj_c, res_c, x); eval_jacobian!(tj_r, res_r, x)
    a_c = @allocated eval_jacobian!(tj_c, res_c, x)
    a_r = @allocated eval_jacobian!(tj_r, res_r, x)
    t_c = timeit(() -> eval_jacobian!(tj_c, res_c, x))
    t_r = timeit(() -> eval_jacobian!(tj_r, res_r, x))
    eval_jacobian!(tj_a, res_a, x)
    t_a = timeit(() -> eval_jacobian!(tj_a, res_a, x))
    t_f = timeit(() -> ForwardDiff.jacobian!(Jf, g!, y, x, cfg))
    t_p = timeit(() -> f(args...)) # plain f on Float64, for reference: the replay has to do about this much work
    best = min(t_c, t_r)
    @printf("   %-28s %dx%-3d dec %4d | compiled %7.2f | replay %7.2f | AUTO %7.2f (chose %-8s %s the better one) | plain f %6.2f | ForwardDiff! %7.2f us   (max diff %.0e)\n",
        name, size(res_c)..., tj_r.traces[1].stats.guards, 1e6 * t_c, 1e6 * t_r, 1e6 * t_a, tj_a.active, t_a <= 1.25 * best ? "is within 25% of" : "is MORE than 25% slower than", 1e6 * t_p, 1e6 * t_f, max(diff, maximum(abs.(res_a .- res_c))))
    println("      auto: ", tj_a.switch_info)
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

bench("f4 (30 conditionals)", f4, (7.3, 0.4))
bench("A \\ b, LU 6x6", (A, b) -> A \ b, (randm(6), randn(rng, 6)))
bench("det and inv, 6x6", A -> vcat(det(A), vec(inv(A))), (randm(6),))
bench("cholesky factor, 8x8", A -> vec(Matrix(cholesky(Symmetric(A' * A + I)).L)), (randn(rng, 8, 8),))
bench("Jacobi SVD 4x4", jacobi_svdvals, (randm(4),))
bench("QuadGK cos(a x) on [0,5]", a -> quadgk(x -> cos(a * x), 0.0, 5.0)[1], (20.0,))
bench("StaticArrays eigvals sym 3x3", a -> collect(eigvals(Symmetric(SMatrix{3,3}(a...)))), (vec(Matrix(Symmetric(randn(rng, 3, 3)))),))

# ============================================================ 3. a cache with several traces
println("\n3. a cache with several traces: average time per call when the calls cycle through inputs that hit different traces")
println("   (compiled guards try the cached traces one after another; replay runs f once and looks the trace up by its decisions)")
function multi(name, f, mk, ninputs)
    inputs = [mk() for _ in 1:ninputs]
    xs = [TP.flatten(a) for a in inputs]
    out = String[]
    for mode in (:compiled, :replay, :auto)
        tj = TracingJacobian(f; in_place=true, guard_mode=mode, max_traces=500)
        res = similar(TracingJacobian(f; guard_mode=:compiled)(inputs[1]...))
        for a in inputs
            tj(res, a...)
        end
        t = timeit(() -> (for x in xs
            eval_jacobian!(tj, res, x)
        end)) / length(xs)
        push!(out, @sprintf("%-8s %6.2f us (%2d traces%s)", string(mode), 1e6 * t, length(tj.traces), mode === :auto ? ", chose " * string(tj.active) : ""))
    end
    println("   ", rpad(name, 32), join(out, " | "))
end
for idx in (5, 6, 7, 8, 9)
    name, f = bank[idx]
    multi(name, f, args_for(name), idx == 6 ? 12 : 40)
end

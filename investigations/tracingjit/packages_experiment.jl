# Try the tracing prototype on functions from Julia numerical packages.
#   julia --project=<env with the packages> investigations/tracingjit/packages_experiment.jl <case name>
# One case per process so that a slow or failing case does not affect the others.

include(joinpath(@__DIR__, "prototype.jl"))
using .TracingPrototype
using ForwardDiff, FiniteDifferences
using Random, Printf, LinearAlgebra
using GenericLinearAlgebra, Roots, QuadGK, LogExpFunctions, StatsFuns, Interpolations, StaticArrays, Rotations, SpecialFunctions

const TP = TracingPrototype
const casename = isempty(ARGS) ? "list" : ARGS[1]
rng = MersenneTwister(7)

# SpecialFunctions: FD has these as Node operations, so a Tracked method is one line each
for f in (:erf, :erfc, :gamma, :loggamma, :digamma)
    @eval SpecialFunctions.$f(a::Tracked) = Tracked(SpecialFunctions.$f(a.val), SpecialFunctions.$f(a.node))
end

# ---------------------------------------------------------------------------- reference derivatives
function make_g(f, args, wrt)
    return function (z)
        pos = 0
        newargs = ntuple(length(args)) do i
            if i in wrt
                len = args[i] isa Number ? 1 : length(args[i])
                piece = z[(pos+1):(pos+len)]
                pos += len
                args[i] isa Number ? piece[1] : reshape(piece, size(args[i]))
            else
                args[i]
            end
        end
        out = f(newargs...)
        return out isa Number ? [out] : collect(vec(out))
    end
end
wrt_vector(args, wrt) = Float64[v for i in wrt for v in (args[i] isa Number ? (args[i],) : vec(args[i]))]

function reference(f, args, wrt)
    g = make_g(f, args, wrt)
    x0 = wrt_vector(args, wrt)
    try
        return ForwardDiff.jacobian(g, x0), "ForwardDiff"
    catch e
        try
            return FiniteDifferences.jacobian(central_fdm(5, 1), g, x0)[1], "finite differences (ForwardDiff failed: " * first(split(sprint(showerror, e), '\n')) * ")"
        catch e2
            return nothing, "none (" * first(split(sprint(showerror, e2), '\n')) * ")"
        end
    end
end

perturbed(args, rel) = map(a -> a isa Number ? a * (1 + rel * randn(rng)) : a .* (1 .+ rel .* randn(rng, size(a)...)), args)
errline(e) = sprint(showerror, e) |> s -> first(split(s, '\n')) |> s -> length(s) > 230 ? s[1:230] * "..." : s

function run_case(name, f, args; wrt=Int[], extra_points=Any[], nperturb=15, reps=2000)
    wrt = isempty(wrt) ? collect(1:length(args)) : collect(wrt)
    println("="^100)
    println("CASE ", name)
    # 1. the plain function
    try
        out = f(args...)
        println("  plain f on Float64: ok (output ", out isa Number ? "scalar" : string(size(out)), ")")
    catch e
        println("  plain f on Float64 FAILED: ", errline(e))
        return
    end
    # TJ_SIZE_ONLY=1: record and differentiate but do not compile, to see how big the trace and its Jacobian are
    if get(ENV, "TJ_SIZE_ONLY", "") == "1"
        tjs = TracingJacobian(f; wrt=wrt)
        TP.init!(tjs, args)
        tr = TP.record!(tjs, args, TP.flatten(args))
        println("  (size only, nothing compiled) ", TP.trace_report(tr))
        return
    end
    # 2. reference derivative
    ref, how = reference(f, args, wrt)
    println("  reference derivative: ", how)
    # 3. trace
    tj = TracingJacobian(f; wrt=wrt, max_traces=200)
    J = nothing
    t_first = @elapsed try
        J = tj(args...)
    catch e
        println("  TRACING FAILED: ", errline(e))
        st = stacktrace(catch_backtrace())
        for fr in st[1:min(6, length(st))]
            println("      at ", fr)
        end
        return
    end
    tr = tj.traces[1]
    println("  first call (record, differentiate, generate, compile): ", round(t_first, digits=2), " s")
    println("  ", TP.trace_report(tr))
    if ref !== nothing
        err = maximum(abs.(J .- ref) ./ max.(abs.(ref), 1.0))
        @printf("  max relative error vs reference: %.2e   (Jacobian is %d x %d)\n", err, size(J)...)
    end
    # 4. nearby inputs: do they reuse the trace?
    h0, m0 = tj.hits, tj.misses
    ok = 0
    for _ in 1:nperturb
        a2 = perturbed(args, 1e-3)
        try
            tj(a2...)
            ok += 1
        catch e
            println("  perturbed call failed: ", errline(e))
            break
        end
    end
    println("  $nperturb perturbed inputs (relative noise 1e-3): ", tj.hits - h0, " cache hits, ", tj.misses - m0, " new traces recorded; ", summary_line(tj))
    # 5. speed of a hit versus ForwardDiff
    try
        tj(args...)
        t_hit = @elapsed for _ in 1:reps
            tj(args...)
        end
        t_hit /= reps
        g = make_g(f, args, wrt)
        x0 = wrt_vector(args, wrt)
        t_fd = Inf
        try
            ForwardDiff.jacobian(g, x0)
            nfd = max(3, min(reps, round(Int, 0.5 / max(1e-7, @elapsed ForwardDiff.jacobian(g, x0)))))
            t_fd = @elapsed(for _ in 1:nfd
                ForwardDiff.jacobian(g, x0)
            end) / nfd
        catch
        end
        @printf("  one call on the cached trace: %.2f us;  ForwardDiff.jacobian: %s\n", 1e6 * t_hit, isfinite(t_fd) ? @sprintf("%.2f us  (cached trace is %.1fx %s)", 1e6 * t_fd, max(t_fd / t_hit, t_hit / t_fd), t_fd > t_hit ? "faster" : "slower") : "not available")
    catch e
        println("  timing failed: ", errline(e))
    end
    # 6. any extra points (other branches)
    for pts in extra_points
        h1, m1 = tj.hits, tj.misses
        try
            Jp = tj(pts...)
            refp, _ = reference(f, pts, wrt)
            e = refp === nothing ? NaN : maximum(abs.(Jp .- refp) ./ max.(abs.(refp), 1.0))
            @printf("  extra point: max relative error %.2e, %s\n", e, tj.misses > m1 ? "new trace" : "cache hit")
        catch ex
            println("  extra point failed: ", errline(ex))
        end
    end
end

# ---------------------------------------------------------------------------- the cases
cases = Dict{String,Function}()
randm(n) = randn(rng, n, n) + n * I
spd(n) = (A = randn(rng, n, n); A' * A + n * I)

# ---- SVD by one-sided Jacobi rotations, written the way real codes are: it stops when no rotation is needed and
# skips pairs of columns that are already orthogonal. The trip count and the skipped rotations depend on the data.
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
for n in (3, 4, 5, 6, 8, 10)
    cases["svd_jacobi_$n"] = () -> run_case("one-sided Jacobi SVD, convergence driven, $(n)x$(n): singular values wrt the matrix", jacobi_svdvals, (randm(n),))
end

# scalar function of the SVD: the nuclear norm (sum of singular values). Its gradient has only n^2 entries, not n x n^2.
for n in (6, 10)
    cases["nuclear_$n"] = () -> run_case("nuclear norm (sum of singular values) from the Jacobi SVD, $(n)x$(n): gradient wrt the matrix", A -> sum(jacobi_svdvals(A)), (randm(n),))
end

# ---- the standard library's generic (non LAPACK) algorithms
cases["lu_solve_6"] = () -> run_case("A \\ b with LinearAlgebra's generic LU, 6x6", (A, b) -> A \ b, (randm(6), randn(rng, 6)))
cases["det_inv_6"] = () -> run_case("det(A) and inv(A), 6x6 (generic LU with pivoting)", A -> vcat(det(A), vec(inv(A))), (randm(6),))
cases["cholesky_8"] = () -> run_case("cholesky(A'A + I) factor, 8x8", A -> vec(Matrix(cholesky(Symmetric(A' * A + I)).L)), (randn(rng, 8, 8),))
cases["qr_6"] = () -> run_case("qr(A) generic Householder, abs of diag(R), 6x6", A -> abs.(diag(qr(A).R)), (randm(6),))

# ---- packages
cases["gla_svdvals_6"] = () -> (run_case("GenericLinearAlgebra.svdvals, 6x6", A -> svdvals(A), (randm(6),)))
cases["gla_svdvals_10"] = () -> (run_case("GenericLinearAlgebra.svdvals, 10x10", A -> svdvals(A), (randm(10),)))
cases["gla_eigvals_8"] = () -> (run_case("GenericLinearAlgebra.eigvals of a symmetric 8x8", A -> eigvals(Symmetric(A + A')), (randn(rng, 8, 8),)))
cases["roots_newton"] = () -> (run_case("Roots.Newton: root of x^3 - a x - 1 as a function of a (data dependent iteration count)",
    a -> find_zero((x -> x^3 - a * x - 1, x -> 3x^2 - a), 1.5 + 0 * a, Roots.Newton()), (2.0,); extra_points=Any[(0.5,), (6.0,), (30.0,)]))
cases["quadgk"] = () -> (run_case("QuadGK.quadgk of exp(-a x^2) on [0, 1]: adaptive subdivision", a -> quadgk(x -> exp(-a * x^2), 0.0, 1.0)[1], (3.0,); extra_points=Any[(30.0,), (300.0,)]))
cases["quadgk_osc"] = () -> (run_case("QuadGK.quadgk of cos(a x) on [0, 5]: more subdivisions", a -> quadgk(x -> cos(a * x), 0.0, 5.0)[1], (20.0,); extra_points=Any[(2.0,), (60.0,)]))
cases["logsumexp_8"] = () -> (run_case("LogExpFunctions.logsumexp, 8 values", x -> logsumexp(x), (randn(rng, 8),)))
cases["log1pexp"] = () -> (run_case("LogExpFunctions.log1pexp (a different formula in each range of x)", x -> log1pexp(x), (-40.0,);
    extra_points=Any[(-5.0,), (5.0,), (20.0,), (40.0,)]))
cases["statsfuns"] = () -> (run_case("StatsFuns normlogpdf, normcdf, gammalogpdf",
    (m, s, x, k) -> [normlogpdf(m, s, x), normcdf(m, s, x), gammalogpdf(k, 2.0, x * x)], (0.3, 1.2, 0.9, 2.5)))
cases["interp_cubic"] = () -> (run_case("Interpolations cubic B-spline of fixed data, evaluated at t (t wrt)",
    t -> (itp = interpolate(sin.(0.7 .* (1:12)), BSpline(Cubic(Line(OnGrid())))); itp(t)), (4.3,); extra_points=Any[(1.2,), (7.9,), (11.5,)]))
cases["interp_cubic_y"] = () -> (run_case("Interpolations cubic B-spline: derivative wrt the data y and t",
    (y, t) -> interpolate(y, BSpline(Cubic(Line(OnGrid()))))(t), (sin.(0.7 .* (1:12)), 4.3)))
cases["static_eigvals3"] = () -> (run_case("StaticArrays eigvals of a symmetric 3x3",
    a -> collect(eigvals(Symmetric(SMatrix{3,3}(a...)))), (vec(Matrix(Symmetric(randn(rng, 3, 3)))),)))
cases["rot_quat"] = () -> (run_case("Rotations: QuatRotation(RotXYZ(r)) (branches on the largest quaternion component)",
    r -> (q = QuatRotation(RotXYZ(r[1], r[2], r[3])); collect(Rotations.params(q))), ([0.3, -0.2, 0.5],); extra_points=Any[([2.8, 0.4, -0.3],), ([0.2, 3.0, 0.1],)]))
cases["rot_angleaxis"] = () -> (run_case("Rotations: AngleAxis(RotXYZ(r)) angle and axis",
    r -> (aa = AngleAxis(RotXYZ(r[1], r[2], r[3])); [rotation_angle(aa), aa.axis_x, aa.axis_y, aa.axis_z]), ([0.3, -0.2, 0.5],); extra_points=Any[([1.0e-3, 2.0e-3, 0.0],)]))

if casename == "list"
    println(join(sort(collect(keys(cases))), "\n"))
else
    haskey(cases, casename) || error("unknown case $casename")
    cases[casename]()
end

# Run the tracing prototype on the hard examples from the conditionals investigation and compare with ForwardDiff.
#   julia --project=. investigations/tracingjit/examples.jl      (ForwardDiff must be on the load path)

include(joinpath(@__DIR__, "prototype.jl"))
using .TracingPrototype
using ForwardDiff
using Random
using Printf

# ---- reference: ForwardDiff on the plain Julia function, differentiating with respect to the arguments in `wrt`
function reference_jacobian(f, args, wrt)
    ranges = TracingPrototype.arg_ranges(args)
    xw = Float64[]
    for i in wrt
        append!(xw, args[i] isa Number ? [args[i]] : vec(collect(Float64, args[i])))
    end
    function g(z)
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
        return out isa Number ? [out] : vec(out)
    end
    return ForwardDiff.jacobian(g, xw)
end

function compare!(label, tj, f, argsets, wrt)
    worst = 0.0
    for args in argsets
        got = tj(args...)
        ref = reference_jacobian(f, args, wrt)
        worst = max(worst, maximum(abs.(got .- ref) ./ max.(abs.(ref), 1.0)))
    end
    @printf("  %-46s max relative error %.1e   %s\n", label, worst, summary_line(tj))
    return worst
end

section(s) = println("\n", "="^100, "\n", s, "\n", "="^100)

rng = MersenneTwister(2026)

# ------------------------------------------------------------------------------------------ 1. basic
section("1. f(a, x) = a < 0 ? x^2 : x^3   (a is a runtime input; one trace per branch)")
f1(a, x) = a < 0 ? x^2 : x^3
compare!("f1", TracingJacobian(f1; wrt=(2,)), f1, [(-1.0, 3.0), (2.0, 3.0), (-0.5, 1.2), (5.0, -2.0), (-3.0, 0.7)], (2,))

# ------------------------------------------------------------------------------------------ 2. max/min/ifelse/abs/clamp
section("2. conditionals written as max, min, ifelse, abs, clamp (no special handling is needed)")
bank = [
    ("min(3t-2, t^2)", (a, t) -> min(3t - 2, t^2)),
    ("max(3t-2, t^2)", (a, t) -> max(3t - 2, t^2)),
    ("max(t, a)^2", (a, t) -> max(t, a)^2),
    ("relu: max(t, 0) * a", (a, t) -> max(t, 0) * a),
    ("abs(t)^2", (a, t) -> abs(t)^2),
    ("clamp(t, -a, a) * t", (a, t) -> clamp(t, -a, a) * t),
    ("ifelse(t < a, t^2, 3t)", (a, t) -> ifelse(t < a, t^2, 3t)),
    ("nested && / ||", (a, t) -> (t < 0 && a > 1) ? t^2 : (t > 1 || a < 0) ? exp(t) : t * a),
]
pts = [(rand(rng) * 4 - 2, rand(rng) * 6 - 3) for _ in 1:12]
for (name, f) in bank
    compare!(name, TracingJacobian(f; wrt=(2,)), f, pts, (2,))
end

# ------------------------------------------------------------------------------------------ 3. guard clause that throws
section("3. guard clause that throws.  f(x) = (x < 0 && throw(DomainError(x)); sqrt(x) * x)")
println("  The eager path enumeration of the IR based implementation follows the throwing side while tracing and fails.")
println("  Here only the path that is actually executed is recorded.")
f3(x) = (x < 0 && throw(DomainError(x, "x must be non negative")); sqrt(x) * x)
tj3 = TracingJacobian(f3)
compare!("f3 on x > 0", tj3, f3, [(2.0,), (5.0,), (0.3,)], (1,))
try
    tj3(-1.0)
    println("  f3(-1.0): no error (unexpected)")
catch e
    println("  f3(-1.0) raises the user's own error, as plain Julia does: ", sprint(showerror, e))
end

# ------------------------------------------------------------------------------------------ 4. many conditionals
section("4. 30 conditionals in a loop (2^30 paths).  u = u > 1.5 ? log(u) : u^2 + p, 30 times")
println("  The eager implementation enumerates every path and stops at its 1024 path limit; tracing only records the paths that run.")
function f4(u, p)
    for i in 1:30
        u = u > 1.5 ? log(u) : u^2 + p
    end
    return u
end
tj4 = TracingJacobian(f4; max_traces=24)
pts4 = [(rand(rng) * 40 + 0.1, rand(rng)) for _ in 1:40]
compare!("f4, 40 random inputs (LRU cap of 24 traces)", tj4, f4, pts4, (1, 2))
t_hit = @elapsed for args in pts4[end-5:end]
    tj4(args...)
end
@printf("  the 6 most recent inputs again (mostly cache hits): %.2f ms for the 6 calls\n", 1000 * t_hit)

# ------------------------------------------------------------------------------------------ 5. while loop, unrolled tracing
section("5. while loop with a data dependent trip count: Newton iteration for sqrt(a), traced as a whole function")
function newton(a)
    s = 1.0 + 0 * a
    while abs(s * s - a) > 1e-12
        s = (s + a / s) / 2
    end
    return s
end
tj5 = TracingJacobian(newton)
as = [2.0, 3.0, 10.0, 123.4, 1e4, 2.0, 7.7, 1e6]
compare!("newton, unrolled (one trace per trip count)", tj5, newton, [(a,) for a in as], (1,))
println("  d sqrt(a)/da at a = 2:  traced = ", tj5(2.0)[1], "   exact = ", 1 / (2 * sqrt(2.0)))

# ------------------------------------------------------------------------------------------ 6. loop mode: Jacobian chaining
section("6. loop mode: trace only the loop BODY and multiply the per iteration Jacobians at runtime")
newton_body(s, p) = [(s[1] + p[1] / s[1]) / 2]
newton_cond(s, p) = abs(s[1]^2 - p[1]) > 1e-12
body_tj = TracingJacobian(newton_body)
for a in [2.0, 10.0, 123.4, 1e4, 1e8]
    s, D, iters = loop_jacobian(body_tj, newton_body, newton_cond, [1.0], [a])
    function ref_loop(z)
        s = z[1:1]
        p = z[2:2]
        while newton_cond(s, p)
            s = newton_body(s, p)
        end
        return s
    end
    ref = ForwardDiff.jacobian(ref_loop, [1.0, a])
    @printf("  a = %-8g iterations = %2d   d s/d a = %.12g   ForwardDiff = %.12g   exact = %.12g\n", a, iters, D[1, 2], ref[1, 2], 1 / (2 * sqrt(a)))
end
println("  body traces cached after all of these (the code size does not depend on the trip count): ", length(body_tj.traces))

println("\n  Same idea for the 30 conditional loop, now with 1000 iterations (unrolled tracing would need a 1000 guard trace):")
cbody(s, p) = [s[1] > 1.5 ? log(s[1]) : s[1]^2 + p[1], s[2] + 1]
ccond(s, p) = s[2] < 1000
cbody_tj = TracingJacobian(cbody)
for (u0, p) in [(7.3, 0.4), (0.2, 0.9), (31.0, 0.1)]
    s, D, iters = loop_jacobian(cbody_tj, cbody, ccond, [u0, 0.0], [p])
    function ref_c(z)
        s = z[1:2]
        pp = z[3:3]
        while ccond(s, pp)
            s = cbody(s, pp)
        end
        return s
    end
    ref = ForwardDiff.jacobian(ref_c, [u0, 0.0, p])
    err = maximum(abs.(D[1, [1, 3]] .- ref[1, [1, 3]]) ./ max.(abs.(ref[1, [1, 3]]), 1e-300))
    @printf("  u0 = %-5g p = %-4g  d u/d(u0, p) = (%.6g, %.6g)   relative difference from ForwardDiff = %.1e\n", u0, p, D[1, 1], D[1, 3], err)
end
println("  body traces cached: ", length(cbody_tj.traces), " (one for each branch of the body)")

# ------------------------------------------------------------------------------------------ 7. B-spline
section("7. cubic B-spline with the curve parameter u known only at runtime (index from the integer part of u)")
bspline_weights(t) = ((1 - t)^3 / 6, (3t^3 - 6t^2 + 4) / 6, (-3t^3 + 3t^2 + 3t + 1) / 6, t^3 / 6)
segment(u) = floor(Int, u)
segment(u::ForwardDiff.Dual) = floor(Int, ForwardDiff.value(u))
function bspline(u, P)
    i = segment(u)
    w = bspline_weights(u - i)
    return sum(w[j] * P[i+j] for j in 1:4)
end
P0 = rand(rng, 12)
us = [0.0, 0.25, 0.9999, 1.0, 1.5, 2.75, 4.5, 6.0, 8.5, 8.99, 0.6, 2.9, 8.5]
tj7 = TracingJacobian(bspline; wrt=(2,))
compare!("d/d(control points), 13 values of u", tj7, bspline, [(u, P0) for u in us], (2,))
tj7b = TracingJacobian(bspline; wrt=(1, 2))
compare!("d/d(u, control points)", tj7b, bspline, [(u, P0) for u in us], (1, 2))
println("  distinct segments among the u values: ", length(unique(floor.(Int, us))), "  (each is a trace; this prototype does not yet share code between them)")

section("   and the vector example  f(u, x) = sum(x[i] for i in Int(mod(u,4))*4 : Int(mod(u,4))*4+4)")
fvec(u, x) = sum(x[i] for i in Int(mod(u, 4))*4:Int(mod(u, 4))*4+4)
tjv = TracingJacobian(fvec; wrt=(2,))
x0 = rand(rng, 12)
compare!("fvec for u = 1, 2, 5.0, 2.0", tjv, fvec, [(1, x0), (2, x0), (5.0, x0), (2.0, x0)], (2,))
println("  u = 1 and 5.0 give the same index, u = 2 and 2.0 give the same index: expect 2 traces and 2 hits, got ", length(tjv.traces), " traces and ", tjv.hits, " hits")
# the guard is on mod(u, 4), so a derivative that does depend on the mod value also works
fmod(u, x) = mod(u, 4) * x[1]^2 + x[2]
compare!("fmod (mod value used in the output), u wrt too", TracingJacobian(fmod; wrt=(1, 2)), fmod, [(1.5, [2.0, 3.0]), (5.5, [2.0, 3.0]), (9.5, [1.0, 1.0])], (1, 2))

# ------------------------------------------------------------------------------------------ cost of a hit and of a miss
section("   cost per call (this machine, one run): a cache hit, a recording, and ForwardDiff for comparison")
function per_call(thunk, reps)
    thunk()
    return @elapsed(for _ in 1:reps
        thunk()
    end) / reps
end
tj_hit = TracingJacobian(f4)
tj_hit(7.3, 0.4) # record
t_hit = per_call(() -> tj_hit(7.3, 0.4), 2000)
t_fd = per_call(() -> ForwardDiff.jacobian(z -> [f4(z[1], z[2])], [7.3, 0.4]), 2000)
@printf("  f4 (30 conditionals): cache hit %.1f us    ForwardDiff %.1f us\n", 1e6 * t_hit, 1e6 * t_fd)
tj_miss = TracingJacobian(f4; max_traces=4)
pts_new = [(rand(rng) * 40 + 0.1, rand(rng)) for _ in 1:30]
tj_miss(pts_new[1]...) # warm up
t_miss = @elapsed for a in pts_new[2:end]
    tj_miss(a...)
end
@printf("  f4: average over %d random inputs (%s): %.1f ms per call, i.e. about %.0fx a ForwardDiff call\n", length(pts_new) - 1, summary_line(tj_miss), 1000 * t_miss / (length(pts_new) - 1), t_miss / (length(pts_new) - 1) / t_fd)
tj_f1 = TracingJacobian(f1; wrt=(2,))
tj_f1(-1.0, 3.0); tj_f1(2.0, 3.0)
t1 = per_call(() -> tj_f1(2.0, 3.0), 20000)
t1fd = per_call(() -> ForwardDiff.derivative(t -> f1(2.0, t), 3.0), 20000)
@printf("  f1 (one conditional): cache hit %.1f us    ForwardDiff %.2f us\n", 1e6 * t1, 1e6 * t1fd)

# ------------------------------------------------------------------------------------------ 8. rand
section("8. a branch on rand(): not fixed by tracing")
function f8(u)
    for i in 1:3
        u = rand() > 0.5 ? u^2 : sqrt(u) + 1
    end
    return u
end
tj8 = TracingJacobian(f8)
println("  traced derivative, 4 calls:  ", [round(tj8(2.0)[1], digits=4) for _ in 1:4], "   <- the same path every call")
println("  ForwardDiff, 4 calls:        ", [round(ForwardDiff.derivative(f8, 2.0), digits=4) for _ in 1:4], "   <- fresh random numbers")

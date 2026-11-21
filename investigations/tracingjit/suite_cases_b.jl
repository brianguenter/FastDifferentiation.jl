# Cases, part B: interpolation, quadrature, root finding, iterative solvers, optimization, special functions, distributions.

const Y10 = sin.(0.6 .* (1:10))
const T8 = collect(1.0:8.0)
const XG = [0.0, 0.5, 1.2, 2.0, 3.1, 4.0, 5.5]
const Z5 = [sin(0.7 * i) * cos(0.5 * j) for i in 1:5, j in 1:5]

# ------------------------------------------------------------------------------------------------ Interpolations
const IP = "Interpolations"
bs(m) = Interpolations.BSpline(m)
add!("itp linear: wrt t", IP, t -> Interpolations.interpolate(Y10, bs(Interpolations.Linear()))(t), R -> (2.0 + 6 * rand(R),); extra=[R -> (1.3,), R -> (9.2,)])
add!("itp quadratic: wrt t", IP, t -> Interpolations.interpolate(Y10, bs(Interpolations.Quadratic(Interpolations.Reflect(Interpolations.OnCell()))))(t), R -> (2.0 + 6 * rand(R),); extra=[R -> (1.3,)])
add!("itp cubic spline: wrt y and t", IP, (y, t) -> Interpolations.interpolate(y, bs(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))(t), R -> (Y10 .+ 0.1 * randn(R, 10), 4.3))
add!("itp constant: wrt y", IP, (y, t) -> Interpolations.interpolate(y, bs(Interpolations.Constant()))(t), R -> (Y10 .+ 0.1 * randn(R, 10), 4.3))
add!("itp gridded linear (nonuniform xs)", IP, (y, t) -> Interpolations.interpolate((XG,), y, Interpolations.Gridded(Interpolations.Linear()))(t), R -> (randn(R, 7), 2.4); extra=[R -> (randn(R, 7), 0.2)])
add!("itp linear with flat extrapolation", IP, t -> Interpolations.extrapolate(Interpolations.interpolate(Y10, bs(Interpolations.Linear())), Interpolations.Flat())(t), R -> (4.5,); extra=[R -> (11.5,), R -> (-0.3,)])
add!("itp linear with line extrapolation", IP, t -> Interpolations.extrapolate(Interpolations.interpolate(Y10, bs(Interpolations.Linear())), Interpolations.Line())(t), R -> (4.5,); extra=[R -> (11.5,)])
add!("itp scaled cubic", IP, t -> Interpolations.scale(Interpolations.interpolate(Y10, bs(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid())))), range(0.0, 4.5; length=10))(t), R -> (2.2,))
add!("itp 2D linear: wrt s and t", IP, (s, t) -> Interpolations.interpolate(Z5, bs(Interpolations.Linear()))(s, t), R -> (1.5 + 2 * rand(R), 2.5 + 2 * rand(R)))
add!("itp 2D cubic: wrt s and t", IP, (s, t) -> Interpolations.interpolate(Z5, bs(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))(s, t), R -> (1.5 + 2 * rand(R), 2.5 + 2 * rand(R)))

# ------------------------------------------------------------------------------------------------ DataInterpolations
const DI = "DataInterpolations"
add!("DI LinearInterpolation", DI, (u, t) -> DataInterpolations.LinearInterpolation(u, T8)(t), R -> (randn(R, 8), 3.4); extra=[R -> (randn(R, 8), 6.9)])
add!("DI QuadraticInterpolation", DI, (u, t) -> DataInterpolations.QuadraticInterpolation(u, T8)(t), R -> (randn(R, 8), 3.4))
add!("DI CubicSpline", DI, (u, t) -> DataInterpolations.CubicSpline(u, T8)(t), R -> (randn(R, 8), 3.4))
add!("DI AkimaInterpolation", DI, (u, t) -> DataInterpolations.AkimaInterpolation(u, T8)(t), R -> (randn(R, 8), 3.4))
add!("DI PCHIPInterpolation", DI, (u, t) -> DataInterpolations.PCHIPInterpolation(u, T8)(t), R -> (randn(R, 8), 3.4))
add!("DI QuadraticSpline", DI, (u, t) -> DataInterpolations.QuadraticSpline(u, T8)(t), R -> (randn(R, 8), 3.4))
add!("DI ConstantInterpolation", DI, (u, t) -> DataInterpolations.ConstantInterpolation(u, T8)(t), R -> (randn(R, 8), 3.4))
add!("DI LagrangeInterpolation", DI, (u, t) -> DataInterpolations.LagrangeInterpolation(u, T8)(t), R -> (randn(R, 8), 3.4))

# ------------------------------------------------------------------------------------------------ BSplineKit
const BK = "BSplineKit"
add!("BSplineKit cubic interpolation", BK, (ys, x) -> BSplineKit.interpolate(collect(0.0:7.0), ys, BSplineKit.BSplineOrder(4))(x), R -> (randn(R, 8), 3.3))
add!("BSplineKit spline derivative", BK, (ys, x) -> (BSplineKit.Derivative(1) * BSplineKit.interpolate(collect(0.0:7.0), ys, BSplineKit.BSplineOrder(4)))(x), R -> (randn(R, 8), 3.3))

# ------------------------------------------------------------------------------------------------ QuadGK and HCubature
const QG = "QuadGK"
add!("quadgk exp(-a x^2) on [0,1]", QG, a -> QuadGK.quadgk(x -> exp(-a * x^2), 0.0, 1.0)[1], R -> (3.0,); extra=[R -> (30.0,), R -> (300.0,)])
add!("quadgk cos(a x) on [0,5]", QG, a -> QuadGK.quadgk(x -> cos(a * x), 0.0, 5.0)[1], R -> (20.0,); extra=[R -> (2.0,)])
add!("quadgk Gaussian on (-Inf, Inf)", QG, a -> QuadGK.quadgk(x -> exp(-a * x^2), -Inf, Inf)[1], R -> (1.7,))
add!("quadgk x^(-p) singular at 0", QG, p -> QuadGK.quadgk(x -> x^(-p), 0.0, 1.0)[1], R -> (0.4,))
add!("quadgk |x-a| kink", QG, a -> QuadGK.quadgk(x -> abs(x - a), 0.0, 1.0)[1], R -> (0.37,))
add!("quadgk with a breakpoint", QG, a -> QuadGK.quadgk(x -> sin(a * x), 0.0, 1.0, 2.0)[1], R -> (3.1,))
add!("quadgk vector integrand", QG, a -> QuadGK.quadgk(x -> [a * x, a * x^2], 0.0, 1.0)[1], R -> (1.3,))
add!("quadgk log(x) exp(-a x)", QG, a -> QuadGK.quadgk(x -> log(x) * exp(-a * x), 0.0, 1.0)[1], R -> (1.5,))
add!("quadgk with tolerance 1e-3", QG, a -> QuadGK.quadgk(x -> exp(-a * x^2), 0.0, 1.0; rtol=1e-3)[1], R -> (3.0,))
const HC = "HCubature"
add!("hcubature 2D Gaussian", HC, a -> HCubature.hcubature(x -> exp(-a * (x[1]^2 + x[2]^2)), (0.0, 0.0), (1.0, 1.0))[1], R -> (2.0,))
add!("hcubature 3D polynomial", HC, a -> HCubature.hcubature(x -> a * x[1] * x[2]^2 + x[3], (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))[1], R -> (2.0,))

# ------------------------------------------------------------------------------------------------ Roots
const RT = "Roots"
f3(a) = x -> x^3 - a * x - 1
d3(a) = x -> 3x^2 - a
add!("Roots Newton", RT, a -> Roots.find_zero((f3(a), d3(a)), 1.5 + 0 * a, Roots.Newton()), R -> (2.0,); extra=[R -> (0.5,), R -> (30.0,)])
add!("Roots Halley", RT, a -> Roots.find_zero((f3(a), d3(a), x -> 6x), 1.5 + 0 * a, Roots.Halley()), R -> (2.0,))
add!("Roots Bisection (bracket)", RT, a -> Roots.find_zero(f3(a), (1.0 + 0 * a, 2.0 + 0 * a), Roots.Bisection()), R -> (2.0,))
add!("Roots Brent", RT, a -> Roots.find_zero(f3(a), (1.0 + 0 * a, 2.0 + 0 * a), Roots.Brent()), R -> (2.0,))
add!("Roots A42", RT, a -> Roots.find_zero(f3(a), (1.0 + 0 * a, 2.0 + 0 * a), Roots.A42()), R -> (2.0,))
add!("Roots Order2 (Steffensen)", RT, a -> Roots.find_zero(f3(a), 1.5 + 0 * a, Roots.Order2()), R -> (2.0,))
add!("Roots Order5", RT, a -> Roots.find_zero(f3(a), 1.5 + 0 * a, Roots.Order5()), R -> (2.0,))
add!("Roots Order8", RT, a -> Roots.find_zero(f3(a), 1.5 + 0 * a, Roots.Order8()), R -> (2.0,))
add!("Roots default bracket", RT, a -> Roots.find_zero(f3(a), (1.0 + 0 * a, 2.0 + 0 * a)), R -> (2.0,))
add!("Roots cos(x) = a x", RT, a -> Roots.find_zero(x -> cos(x) - a * x, (0.0 + 0 * a, 1.0 + 0 * a), Roots.Brent()), R -> (0.8,))

# ------------------------------------------------------------------------------------------------ iterative linear solvers
const IS = "IterativeSolvers"
sym(R, n) = (A = randn(R, n, n); A' * A + n * I)
add!("IterativeSolvers cg 6x6", IS, (A, b) -> IterativeSolvers.cg(A, b; reltol=1e-10), R -> (sym(R, 6), randn(R, 6)))
add!("IterativeSolvers gmres 6x6", IS, (A, b) -> IterativeSolvers.gmres(A, b; reltol=1e-10), R -> (dd(R, 6), randn(R, 6)))
add!("IterativeSolvers minres 6x6", IS, (A, b) -> IterativeSolvers.minres(A, b; reltol=1e-10), R -> (sym(R, 6), randn(R, 6)))
add!("IterativeSolvers bicgstabl 6x6", IS, (A, b) -> IterativeSolvers.bicgstabl(A, b, 2; reltol=1e-10), R -> (dd(R, 6), randn(R, 6)))
add!("IterativeSolvers lsmr 8x4", IS, (A, b) -> IterativeSolvers.lsmr(A, b), R -> (randn(R, 8, 4), randn(R, 8)))
add!("IterativeSolvers lsqr 8x4", IS, (A, b) -> IterativeSolvers.lsqr(A, b), R -> (randn(R, 8, 4), randn(R, 8)))
const KR = "Krylov"
add!("Krylov cg 6x6", KR, (A, b) -> Krylov.cg(A, b; rtol=1e-10)[1], R -> (sym(R, 6), randn(R, 6)))
add!("Krylov minres 6x6", KR, (A, b) -> Krylov.minres(A, b; rtol=1e-10)[1], R -> (sym(R, 6), randn(R, 6)))
add!("Krylov gmres 6x6", KR, (A, b) -> Krylov.gmres(A, b; rtol=1e-10)[1], R -> (dd(R, 6), randn(R, 6)))
add!("Krylov bicgstab 6x6", KR, (A, b) -> Krylov.bicgstab(A, b; rtol=1e-10)[1], R -> (dd(R, 6), randn(R, 6)))

# ------------------------------------------------------------------------------------------------ Optim, LsqFit
const OP = "Optim"
rosen(p) = x -> (1 - x[1])^2 + p * (x[2] - x[1]^2)^2
rosen_g!(p) = (g, x) -> (g[1] = -2 * (1 - x[1]) - 4p * x[1] * (x[2] - x[1]^2); g[2] = 2p * (x[2] - x[1]^2); g)
add!("Optim NelderMead (Rosenbrock p)", OP, p -> Optim.minimizer(Optim.optimize(rosen(p), [0.0 + 0 * p, 0.0 + 0 * p], Optim.NelderMead(); iterations=40)), R -> (3.0,))
add!("Optim GradientDescent (analytic gradient)", OP, p -> Optim.minimizer(Optim.optimize(rosen(p), rosen_g!(p), [0.0 + 0 * p, 0.0 + 0 * p], Optim.GradientDescent(); iterations=20)), R -> (3.0,))
add!("Optim BFGS (analytic gradient)", OP, p -> Optim.minimizer(Optim.optimize(rosen(p), rosen_g!(p), [0.0 + 0 * p, 0.0 + 0 * p], Optim.BFGS(); iterations=20)), R -> (3.0,))
add!("Optim quadratic LBFGS", OP, c -> Optim.minimizer(Optim.optimize(x -> sum((x .- c) .^ 2) + 0.1 * sum(x .^ 4), (g, x) -> (g .= 2 .* (x .- c) .+ 0.4 .* x .^ 3), zeros(typeof(c[1]), 3), Optim.LBFGS(); iterations=15)), R -> (randn(R, 3),))
const LF = "LsqFit"
add!("LsqFit exponential fit", LF, ys -> LsqFit.curve_fit((x, p) -> p[1] .* exp.(-p[2] .* x), collect(0.0:0.5:3.0), ys, [1.0, 0.5]).param, R -> (2.0 .* exp.(-0.7 .* collect(0.0:0.5:3.0)) .+ 0.01 * randn(R, 7),))

# ------------------------------------------------------------------------------------------------ differential equations
const OD = "OrdinaryDiffEqTsit5"
add!("Tsit5 exponential decay", OD, (u0, p) -> OrdinaryDiffEqTsit5.solve(OrdinaryDiffEqTsit5.ODEProblem((u, q, t) -> -q[1] * u, u0, (0.0, 1.0), [p]), OrdinaryDiffEqTsit5.Tsit5(); reltol=1e-8, abstol=1e-10, save_everystep=false).u[end], R -> (1.3, 0.8))
add!("Tsit5 Lotka-Volterra", OD, (u0, p) -> OrdinaryDiffEqTsit5.solve(OrdinaryDiffEqTsit5.ODEProblem((u, q, t) -> [q[1] * u[1] - q[2] * u[1] * u[2], -q[3] * u[2] + q[4] * u[1] * u[2]], u0, (0.0, 1.0), p), OrdinaryDiffEqTsit5.Tsit5(); reltol=1e-6, abstol=1e-8, save_everystep=false).u[end], R -> ([1.0, 1.0], [1.5, 1.0, 3.0, 1.0]))
const TI = "TaylorIntegration"
add!("TaylorIntegration decay", TI, (x0, p) -> TaylorIntegration.taylorinteg((x, q, t) -> -q * x, x0, 0.0, 1.0, 15, 1e-16, p)[2][end], R -> (1.3, 0.8))

# ------------------------------------------------------------------------------------------------ special functions
const SP = "SpecialFunctions"
add!("erf and erfc", SP, x -> [SpecialFunctions.erf(x), SpecialFunctions.erfc(x)], R -> (0.4,))
add!("gamma and loggamma", SP, x -> [SpecialFunctions.gamma(x), SpecialFunctions.loggamma(x)], R -> (2.3,))
add!("digamma", SP, x -> SpecialFunctions.digamma(x), R -> (2.3,))
add!("beta and logbeta", SP, (a, b) -> [SpecialFunctions.beta(a, b), SpecialFunctions.logbeta(a, b)], R -> (2.3, 1.7))
add!("erfinv", SP, x -> SpecialFunctions.erfinv(x), R -> (0.3,); expect_failure=true)
add!("trigamma", SP, x -> SpecialFunctions.trigamma(x), R -> (2.3,); expect_failure=true)
add!("gamma_inc", SP, (a, x) -> SpecialFunctions.gamma_inc(a, x)[1], R -> (2.3, 1.1); expect_failure=true)
add!("besselj (Float64 only)", SP, x -> SpecialFunctions.besselj(0.5, x), R -> (2.3,); expect_failure=true)
const BS = "Bessels"
for (nm, fn, xs) in (("besselj0", Bessels.besselj0, (0.7, 4.0, 30.0)), ("besselj1", Bessels.besselj1, (0.7, 4.0, 30.0)), ("bessely0", Bessels.bessely0, (0.7, 4.0, 30.0)),
    ("bessely1", Bessels.bessely1, (0.7, 4.0, 30.0)), ("besseli0", Bessels.besseli0, (0.7, 4.0, 30.0)), ("besselk0", Bessels.besselk0, (0.7, 4.0, 30.0)))
    add!("Bessels $nm (3 ranges)", BS, fn, R -> (xs[1],); extra=[R -> (xs[2],), R -> (xs[3],)])
end
add!("Bessels besselj order 2.5", BS, x -> Bessels.besselj(2.5, x), R -> (3.3,))
add!("Bessels airyai", BS, x -> Bessels.airyai(x), R -> (0.8,); extra=[R -> (-3.0,), R -> (6.0,)])
const LE = "LogExpFunctions"
add!("logsumexp 8", LE, x -> LogExpFunctions.logsumexp(x), R -> (randn(R, 8),))
add!("softmax 6", LE, x -> LogExpFunctions.softmax(x), R -> (randn(R, 6),))
add!("log1pexp (4 ranges)", LE, x -> LogExpFunctions.log1pexp(x), R -> (-40.0,); extra=[R -> (-5.0,), R -> (5.0,), R -> (20.0,), R -> (40.0,)])
add!("log1mexp, logexpm1", LE, x -> [LogExpFunctions.log1mexp(-x), LogExpFunctions.logexpm1(x)], R -> (0.4,); extra=[R -> (1.5,)])
add!("logaddexp, logsubexp", LE, (a, b) -> [LogExpFunctions.logaddexp(a, b), LogExpFunctions.logsubexp(max(a, b), min(a, b))], R -> (1.2, 0.3))
add!("logistic and logit", LE, x -> [LogExpFunctions.logistic(x), LogExpFunctions.logit(0.3 + 0.4 * LogExpFunctions.logistic(x))], R -> (0.4,); extra=[R -> (40.0,), R -> (-40.0,)])
add!("xlogx, xlogy, xlog1py", LE, (x, y) -> [LogExpFunctions.xlogx(x), LogExpFunctions.xlogy(x, y), LogExpFunctions.xlog1py(x, y)], R -> (0.7, 1.9))
add!("cloglog, log1psq, log1pmx", LE, x -> [LogExpFunctions.cloglog(0.3 + 0.4 * LogExpFunctions.logistic(x)), LogExpFunctions.log1psq(x), LogExpFunctions.log1pmx(x / 3)], R -> (0.4,))
const SF = "StatsFuns"
add!("StatsFuns norm pdf cdf", SF, (m, s, x) -> [StatsFuns.normlogpdf(m, s, x), StatsFuns.normcdf(m, s, x), StatsFuns.normlogcdf(m, s, x)], R -> (0.3, 1.2, 0.9))
add!("StatsFuns norminvcdf", SF, p -> StatsFuns.norminvcdf(p), R -> (0.3,); expect_failure=true)
add!("StatsFuns gamma beta logpdf", SF, (k, th, a, b, x) -> [StatsFuns.gammalogpdf(k, th, x), StatsFuns.betalogpdf(a, b, 0.3 + 0.4 * x)], R -> (2.5, 1.2, 2.0, 3.0, 0.5))
add!("StatsFuns t chisq F logpdf", SF, (nu, x) -> [StatsFuns.tdistlogpdf(nu, x), StatsFuns.chisqlogpdf(nu, x * x + 0.1), StatsFuns.fdistlogpdf(nu, nu + 1, x * x + 0.1)], R -> (4.5, 0.8))
add!("StatsFuns poisson binomial", SF, (lam, p) -> [StatsFuns.poislogpmf(lam, 3), StatsFuns.binomlogpdf(10, p, 3)], R -> (2.3, 0.4))

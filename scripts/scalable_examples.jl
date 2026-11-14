# Scalable example expressions shared by scripts/measure_dominators.jl and
# scripts/measure_segments.jl. Each function takes a size and returns the
# vector of output expressions (the roots of the derivative graph).
#
# spherical harmonics, Chebyshev polynomials and the heat equation of issue #60
# come from test/runtests.jl; the Rosenbrock function, the sin*exp chain and the
# dense map are simple scalable functions; the linear MPC rollout is from issue
# #26. `grid` is a k x k grid, used to show a graph with exponentially many paths.

import FastDifferentiation as FD

# The tests use Memoize.jl, which is only a test dependency; memoize by hand here.
const MEMO = IdDict{Any,Any}()  # keyed by identity: hashing a Node (a Real) by value fails
macro memo(def)
    call, body = def.args[1], def.args[2]
    name, args = call.args[1], call.args[2:end]
    argnames = [a isa Expr ? a.args[1] : a for a in args]
    esc(:($call = get!(() -> $body, MEMO, ($(QuoteNode(name)), $(argnames...)))))
end

# spherical harmonics, copied from test/runtests.jl
@memo function P(l, m, z)
    if l == 0 && m == 0
        return 1.0
    elseif l == m
        return (1 - 2m) * P(m - 1, m - 1, z)
    elseif l == m + 1
        return (2m + 1) * z * P(m, m, z)
    else
        return ((2l - 1) / (l - m) * z * P(l - 1, m, z) - (l + m - 1) / (l - m) * P(l - 2, m, z))
    end
end
@memo S(m, x, y) = m == 0 ? 0 : x * C(m - 1, x, y) - y * S(m - 1, x, y)
@memo C(m, x, y) = m == 0 ? 1 : x * S(m - 1, x, y) + y * C(m - 1, x, y)
factorial_approximation(n1) = sqrt(2 * π * n1) * (n1 / ℯ * sqrt(n1 * sinh(1 / n1) + 1 / (810 * n1^6)))^n1
@memo N(l, m) = m == 0 ? sqrt((2l + 1 / (4π))) : sqrt((2l + 1) / 2π * factorial_approximation(l - m) / factorial_approximation(l + m))
@memo Y(l, m, x, y, z) = m < 0 ? N(l, abs(m)) * P(l, abs(m), z) * S(abs(m), x, y) : N(l, m) * P(l, m, z) * C(m, x, y)
function spherical_harmonics(max_l)
    FD.@variables x y z
    return FD.Node[Y(l, m, x, y, z) for l in 0:max_l-1 for m in -l:l]  # as in the tests; Y(0, 0) is a plain number
end

# Chebyshev polynomial, copied from test/runtests.jl
@memo Chebyshev(n, x) = n == 0 ? 1 : n == 1 ? x : 2 * x * Chebyshev(n - 1, x) - Chebyshev(n - 2, x)
function chebyshev(n)
    FD.@variables x
    return [Chebyshev(n, x)]
end

# heat equation, issue #60 (test/runtests.jl, "jacobian_times_v")
function heat(Ns)
    u = FD.make_variables(:u, Ns)
    p = FD.make_variables(:p, 2)
    dx = 1.0 / (Ns - 1)
    return p[1] * (circshift(u, 1) .- 2 * u .+ circshift(u, -1)) / dx^2
end

# simple scalable functions
function rosenbrock(n)
    x = FD.make_variables(:x, n)
    return [sum(100 * (x[i+1] - x[i]^2)^2 + (1 - x[i])^2 for i in 1:n-1)]
end
function sin_exp_chain(n)
    y = FD.make_variables(:y, n)
    return [sum(sin(y[i] * y[i+1]) * exp(y[i]) for i in 1:n-1)]
end
function dense_map(n)
    z = FD.make_variables(:z, n)
    return [z[i] * sin(sum(z[j] * (i + j) / n for j in 1:n)) for i in 1:n]
end

# linear MPC dynamics constraints, issue #26, with a fixed random system (nx = 4, nu = 2)
const MPC_A = [0.31 -0.12 0.05 0.22; 0.08 0.17 -0.26 0.04; -0.19 0.11 0.09 -0.07; 0.02 -0.15 0.28 0.13]
const MPC_B = [0.9 -0.4; 0.3 1.1; -0.7 0.2; 0.5 -0.8]
function mpc(N)
    nx, nu = size(MPC_B)
    x = FD.make_variables(:x, nx, N + 1)
    u = FD.make_variables(:u, nu, N)
    x0 = FD.make_variables(:x0, nx)
    xp = Matrix{FD.Node}(undef, nx, N + 1)
    xp[:, 1] .= x0
    for i in 1:N
        xp[:, i+1] = MPC_A * xp[:, i] + MPC_B * u[:, i]
    end
    return vec(xp .- x)
end

# k x k grid: g[i,j] = g[i-1,j] * g[i,j-1], with sin and cos along the two edges
# of the grid. There are binomial(2k-2, k-1) paths from g[k,k] down to x.
function grid(k)
    FD.@variables x
    g = Matrix{FD.Node}(undef, k, k)
    g[1, 1] = x
    for i in 2:k
        g[i, 1] = sin(g[i-1, 1])
    end
    for j in 2:k
        g[1, j] = cos(g[1, j-1])
    end
    for i in 2:k, j in 2:k
        g[i, j] = g[i-1, j] * g[i, j-1]
    end
    return [g[k, k]]
end

#!/usr/bin/env julia
# Reproduces the dominator-table charts in HOW_FASTDIFFERENTIATION_WORKS.md,
# section 14 (Cost).
#
# Usage (from the repository root):
#     julia --project=. scripts/measure_dominators.jl          # writes img/fig-dom-*.svg
#     julia --project=. scripts/measure_dominators.jl outdir   # writes the charts elsewhere
#
# For each scalable example it builds the derivative graph at several sizes,
# runs the dominator and postdominator tables exactly as
# compute_factorable_subgraphs does (one per root, one per variable), and counts
# the steps taken by the `intersect` walks. The counter is added by redefining
# FastDifferentiation.intersection with the same loop plus a counter.
#
# It writes one chart per example (intersect steps against nodes + edges) and
# prints the measurements as a Markdown table.
#
# Examples: spherical harmonics, Chebyshev polynomials and the heat equation of
# issue #60 (all from test/runtests.jl); the Rosenbrock function, a chain of
# sin(x_i x_{i+1}) exp(x_i) terms and a dense map (from scripts/measure_groups.jl);
# and the linear model-predictive-control rollout of issue #26.

import FastDifferentiation as FD

# The tests use Memoize.jl, which is only a test dependency; memoize by hand here.
const MEMO = IdDict{Any,Any}()  # keyed by identity: hashing a Node (a Real) by value fails
macro memo(def)
    call, body = def.args[1], def.args[2]
    name, args = call.args[1], call.args[2:end]
    argnames = [a isa Expr ? a.args[1] : a for a in args]
    esc(:($call = get!(() -> $body, MEMO, ($(QuoteNode(name)), $(argnames...)))))
end

const STEPS = Ref(0)
@eval FD function intersection(order_test, node1::T, node2::T, idoms::Dict{T,T}, start_index::T) where {T<:Integer}
    while node1 != node2
        $STEPS[] += 1
        if order_test(node1, node2)
            node1 = idoms[node1]
        else
            node2 = idoms[node2]
        end
    end
    return node1
end

# ---------------------------------------------------------------- examples
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

# from scripts/measure_groups.jl
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

const EXAMPLES = [
    ("spherical harmonics", "fig-dom-spherical-harmonics.svg", s -> "order $s", spherical_harmonics, 4:4:24),
    ("Chebyshev polynomial", "fig-dom-chebyshev.svg", s -> "degree $s", chebyshev, 50:50:300),
    ("heat equation (issue #60)", "fig-dom-heat.svg", s -> "$s points", heat, 40:40:240),
    ("Rosenbrock function", "fig-dom-rosenbrock.svg", s -> "$s variables", rosenbrock, 20:20:120),
    ("Σ sin(xᵢxᵢ₊₁)·exp(xᵢ)", "fig-dom-sin-exp-chain.svg", s -> "$s variables", sin_exp_chain, 20:20:120),
    ("dense map xᵢ·sin(Σⱼ xⱼ(i+j)/n)", "fig-dom-dense-map.svg", s -> "n = $s", dense_map, 5:5:30),
    ("linear MPC rollout (issue #26)", "fig-dom-mpc.svg", s -> "horizon $s", mpc, 5:7:40),
]

# ---------------------------------------------------------------- measurement
function measure(roots)
    gr = FD.DerivativeGraph(roots)
    E = sum(length(FD.children(e)) for e in values(FD.edges(gr)))
    V = length(FD.nodes(gr))
    tables = FD.codomain_dimension(gr) + FD.domain_dimension(gr)
    STEPS[] = 0
    doms = Dict{Int,Int}()
    for i in 1:FD.codomain_dimension(gr)
        FD.compute_dom_table(gr, true, i, FD.root_index_to_postorder_number(gr, i), doms)
    end
    for j in 1:FD.domain_dimension(gr)
        FD.compute_dom_table(gr, false, j, FD.variable_index_to_postorder_number(gr, j), doms)
    end
    return (V=V, E=E, tables=tables, steps=STEPS[])
end

# ---------------------------------------------------------------- charts
# Light palette, matching the other figures: the series uses categorical slot 1
# of the dataviz reference palette; text and axes use neutral text tokens.
const SURFACE = "#ffffff"
const SERIES = "#2a78d6"
const TEXT = "#16202b"
const MUTED = "#586473"
const GRID = "#e3e7ec"
const REF = "#9aa4b0"
const SANS = "'Space Grotesk', 'Segoe UI', Helvetica, Arial, sans-serif"
const MONO = "'JetBrains Mono', Consolas, 'DejaVu Sans Mono', Menlo, monospace"

"About five round tick values from 0 to at least `hi`."
function ticks(hi)
    hi <= 0 && return [0]
    raw = hi / 4
    mag = 10.0^floor(log10(raw))
    step = first(s for s in (1, 2, 2.5, 5, 10) .* mag if s >= raw)
    return collect(0:step:(ceil(hi / step) * step))
end
function fmt(v)
    v < 10_000 && return replace(string(round(Int, v)), r"(?<=\d)(?=(\d{3})+$)" => ",")
    k = round(v / 1000, digits=v >= 100_000 ? 0 : 1)
    return (isinteger(k) ? string(Int(k)) : string(k)) * "K"
end
r1(v) = string(round(v, digits=1))

function chart(title, sizephrase, sizes, data)
    W, H = 480, 300
    L, R, T, B = 66, 22, 62, 46
    xs = [d.V + d.E for d in data]; ys = [d.steps for d in data]
    xt = ticks(maximum(xs))
    yt = ticks(1.15 * maximum(ys))   # headroom so the end label clears the legend
    ytop = maximum(ys) > 0 ? last(yt) : 1
    X(v) = L + (W - L - R) * v / last(xt)
    Yc(v) = H - B - (H - T - B) * v / ytop
    p = String[]
    push!(p, """<rect width="$W" height="$H" rx="12" fill="$SURFACE"/>""")
    push!(p, """<text x="16" y="26" font-family="$SANS" font-size="15" font-weight="700" fill="$TEXT">$title</text>""")
    # legend: measured series, and the proportional reference
    push!(p, """<circle cx="22" cy="44" r="4" fill="$SERIES" stroke="$SURFACE" stroke-width="2"/>""")
    push!(p, """<text x="32" y="48" font-family="$SANS" font-size="12" fill="$MUTED">intersect steps</text>""")
    if last(ys) > 0
        push!(p, """<line x1="140" y1="44" x2="162" y2="44" stroke="$REF" stroke-width="1.5" stroke-dasharray="4 3"/>""")
        push!(p, """<text x="168" y="48" font-family="$SANS" font-size="12" fill="$MUTED">proportional to size, through the largest example</text>""")
    end
    for v in yt
        push!(p, """<line x1="$L" y1="$(r1(Yc(v)))" x2="$(W - R)" y2="$(r1(Yc(v)))" stroke="$GRID" stroke-width="1"/>""")
        push!(p, """<text x="$(L - 8)" y="$(r1(Yc(v) + 4))" font-family="$MONO" font-size="11" text-anchor="end" fill="$MUTED">$(fmt(v))</text>""")
    end
    for v in xt
        push!(p, """<line x1="$(r1(X(v)))" y1="$(H - B)" x2="$(r1(X(v)))" y2="$(H - B + 4)" stroke="$MUTED" stroke-width="1"/>""")
        push!(p, """<text x="$(r1(X(v)))" y="$(H - B + 17)" font-family="$MONO" font-size="11" text-anchor="middle" fill="$MUTED">$(fmt(v))</text>""")
    end
    push!(p, """<line x1="$L" y1="$(H - B)" x2="$(W - R)" y2="$(H - B)" stroke="$MUTED" stroke-width="1"/>""")
    push!(p, """<text x="$(L + (W - L - R) / 2)" y="$(H - 8)" font-family="$SANS" font-size="12" text-anchor="middle" fill="$MUTED">size of the derivative graph (nodes + edges)</text>""")
    if last(ys) > 0
        push!(p, """<line x1="$(X(0))" y1="$(Yc(0))" x2="$(r1(X(last(xs))))" y2="$(r1(Yc(last(ys))))" stroke="$REF" stroke-width="1.5" stroke-dasharray="4 3"/>""")
    else
        push!(p, """<text x="$(L + (W - L - R) / 2)" y="$(r1(Yc(0.5)))" font-family="$SANS" font-size="12" text-anchor="middle" fill="$MUTED">no intersect steps at any size</text>""")
    end
    pts = join(["$(r1(X(x))),$(r1(Yc(y)))" for (x, y) in zip(xs, ys)], " ")
    push!(p, """<polyline points="$pts" fill="none" stroke="$SERIES" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>""")
    for (x, y) in zip(xs, ys)
        push!(p, """<circle cx="$(r1(X(x)))" cy="$(r1(Yc(y)))" r="4" fill="$SERIES" stroke="$SURFACE" stroke-width="2"/>""")
    end
    # direct label on the largest example only
    push!(p, """<text x="$(r1(X(last(xs)) - 8))" y="$(r1(Yc(last(ys)) - 10))" font-family="$MONO" font-size="11" text-anchor="end" fill="$TEXT">$(fmt(last(ys))) steps at $(sizephrase(last(sizes)))</text>""")
    return """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $W $H" width="$W" height="$H" role="img" aria-label="$title: intersect steps against graph size">
  $(join(p, "\n  "))
</svg>
"""
end

outdir = isempty(ARGS) ? joinpath(@__DIR__, "..", "img") : ARGS[1]
mkpath(outdir)
measure(chebyshev(5))   # compile before measuring

println("| example | size | nodes + edges | tables (roots + variables) | intersect steps | steps per (node + edge) per table |")
println("|:--|--:|--:|--:|--:|--:|")
for (title, file, sizephrase, f, sizes) in EXAMPLES
    data = [measure(f(s)) for s in sizes]
    for (s, d) in zip(sizes, data)
        println("| $title | $(sizephrase(s)) | $(d.V + d.E) | $(d.tables) | $(d.steps) | $(round(d.steps / ((d.V + d.E) * d.tables), digits=3)) |")
    end
    write(joinpath(outdir, file), chart(title, sizephrase, sizes, data))
end

#!/usr/bin/env julia
# Reproduces the dominator-table measurements in HOW_FASTDIFFERENTIATION_WORKS.md,
# section 14 (Cost).
#
# Usage (from the repository root):
#     julia --project=. scripts/measure_dominators.jl
#
# For the scalable examples used in the tests (spherical harmonics, Chebyshev
# polynomials, the heat equation of issue #60) it builds the derivative graph,
# runs the dominator and postdominator tables exactly as
# compute_factorable_subgraphs does (one per root, one per variable), and counts
# the steps taken by the `intersect` walks. The counter is added by redefining
# FastDifferentiation.intersection with the same loop plus a counter.
#
# `steps per edge per table` is total steps / (edges x tables). If the
# dominator stage were linear in the size of the graph this would stay
# constant as the examples grow.

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

# ---- examples, copied from test/runtests.jl
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

@memo Chebyshev(n, x) = n == 0 ? 1 : n == 1 ? x : 2 * x * Chebyshev(n - 1, x) - Chebyshev(n - 2, x)
function chebyshev(n)
    FD.@variables x
    return [Chebyshev(n, x)]
end

function heat(Ns)
    u = FD.make_variables(:u, Ns)
    p = FD.make_variables(:p, 2)
    dx = 1.0 / (Ns - 1)
    return p[1] * (circshift(u, 1) .- 2 * u .+ circshift(u, -1)) / dx^2
end

# ---- measurement
function measure(label, roots; quiet=false)
    gr = FD.DerivativeGraph(roots)
    E = sum(length(FD.children(e)) for e in values(FD.edges(gr)))
    tables = FD.codomain_dimension(gr) + FD.domain_dimension(gr)
    STEPS[] = 0
    t = @elapsed begin
        doms = Dict{Int,Int}()
        for i in 1:FD.codomain_dimension(gr)
            FD.compute_dom_table(gr, true, i, FD.root_index_to_postorder_number(gr, i), doms)
        end
        for j in 1:FD.domain_dimension(gr)
            FD.compute_dom_table(gr, false, j, FD.variable_index_to_postorder_number(gr, j), doms)
        end
    end
    quiet && return
    println("| ", label, " | ", length(FD.nodes(gr)), " | ", E, " | ", tables, " | ", STEPS[], " | ",
        round(STEPS[] / (E * tables), digits=2), " | ", round(1000t, digits=1), " |")
end

measure("warm-up", chebyshev(5); quiet=true)  # compile before timing
println("| example | nodes | edges | tables (roots + variables) | intersect steps | steps per edge per table | time (ms) |")
println("|:--|--:|--:|--:|--:|--:|--:|")
for l in (5, 10, 15, 20)
    measure("spherical harmonics, order $l", spherical_harmonics(l))
end
for n in (20, 40, 80, 160)
    measure("Chebyshev polynomial, degree $n", chebyshev(n))
end
for Ns in (5, 20, 80, 320)
    measure("heat equation, $Ns points", heat(Ns))
end

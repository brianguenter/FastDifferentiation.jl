# Is my traced graph a worse input to FD than a graph built directly from FD Nodes?
# The same branch-free algorithm (fixed sweeps, atan rotation angle) is evaluated two ways:
#   (A) directly on a matrix of FD variables (what a normal FD user would do)
#   (B) on Tracked numbers (what the tracing prototype does)
# and for each the graph is differentiated with FD.jacobian and FD.reverse_AD.
#   julia --project=<env> investigations/tracingjit/graph_check.jl 3 4 5

include(joinpath(@__DIR__, "prototype.jl"))
using .TracingPrototype
using FastDifferentiation
using Random, Printf, LinearAlgebra

const TP = TracingPrototype
const FD = FastDifferentiation
const Node = FD.Node

function jacobi_fixed(A, sweeps)
    m, n = size(A)
    B = copy(A)
    for sweep in 1:sweeps
        for p in 1:n-1, q in p+1:n
            α = zero(eltype(B)); β = zero(eltype(B)); γ = zero(eltype(B))
            for i in 1:m
                α += B[i, p] * B[i, p]
                β += B[i, q] * B[i, q]
                γ += B[i, p] * B[i, q]
            end
            θ = atan(2γ, β - α) / 2
            c = cos(θ); s = sin(θ)
            for i in 1:m
                bp = B[i, p]; bq = B[i, q]
                B[i, p] = c * bp - s * bq
                B[i, q] = s * bp + c * bq
            end
        end
    end
    return sum(sqrt(sum(abs2, B[:, j])) for j in 1:n)
end

ops(x) = count(FD.is_tree, FD.all_nodes(x))

# number of structurally distinct nodes: nodes with the same operation and structurally identical children are merged.
# If this equals the number of nodes, the graph already shares every repeated subexpression.
function distinct_structures(roots)
    ids = IdDict{Node,Int}()
    table = Dict{Any,Int}()
    function id(n)
        get!(ids, n) do
            key = FD.is_tree(n) ? (FD.value(n), map(id, FD.children(n))...) : (FD.is_variable(n) ? (:var, objectid(n)) : (:const, FD.value(n)))
            get!(table, key, length(table) + 1)
        end
    end
    foreach(id, roots)
    return length(table)
end

function report(label, root, vars)
    t_fd = @elapsed J = FD.jacobian([root], vars)
    t_rev = @elapsed R = FD.reverse_AD(root, vars)
    base = ops([root])
    @printf("   %-34s graph ops %7d (distinct structures %7d, nodes %7d)   FD.jacobian ops %8d (%5.1fx)   reverse_AD ops %7d (%4.1fx)\n",
        label, base, distinct_structures([root]), length(FD.all_nodes([root])), ops(J), ops(J) / base, ops(R), ops(R) / base)
end

rng = MersenneTwister(3)
sizes = isempty(ARGS) ? [3, 4, 5] : parse.(Int, ARGS)
sweeps = 3
for n in sizes
    println("n = $n, $sweeps sweeps")
    A0 = randn(rng, n, n) + n * I

    # (A) built directly from FD variables
    varsA = [Node(Symbol("a", i)) for i in 1:length(A0)]
    rootA = jacobi_fixed(reshape(varsA, n, n), sweeps)
    report("(A) built directly from FD Nodes", rootA, varsA)

    # (B) built by running the same code on Tracked numbers
    varsB = [Node(Symbol("b", i)) for i in 1:length(A0)]
    tracked = reshape([TP.Tracked(A0[i], varsB[i]) for i in 1:length(A0)], n, n)
    rec = TP.Recording()
    TP.CURRENT[] = rec
    rootB = try
        jacobi_fixed(tracked, sweeps).node
    finally
        TP.CURRENT[] = nothing
    end
    report("(B) recorded with Tracked", rootB, varsB)
end

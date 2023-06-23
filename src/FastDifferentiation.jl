module FastDifferentiation

using TermInterface
import SymbolicUtils
using SymbolicUtils: arguments
using StaticArrays
using SpecialFunctions
using NaNMath
# using Statistics
using RuntimeGeneratedFunctions
import Base: iterate
using UUIDs
using SparseArrays
using DataStructures
module AutomaticDifferentiation
struct NoDeriv
end
export NoDeriv
end #module

const INVARIANTS = true

"""This macro is used to create invariant test code that is dependent on the global constant `INVARIANTS`. If `INVARIANTS` is false then the test code will not be inserted into the program and there will be no run time overhead. If `INVARIANTS` is true then the code will be inserted. Code that tests invariants tends to increase run time substantially so only set `INVARIANTS` true when you are debugging or testing."""
macro invariant(ex, msgs...)
    if INVARIANTS
        return :(@assert $(esc(ex)) $(esc(msgs)))
    end
end



RuntimeGeneratedFunctions.init(@__MODULE__)

const DefaultNodeIndexType = Int64

include("Utilities.jl")
include("BitVectorFunctions.jl")
include("ExpressionGraph.jl")
include("PathEdge.jl")
include("DerivativeGraph.jl")
include("Reverse.jl")
include("GraphProcessing.jl")
include("FactorableSubgraph.jl")
include("Factoring.jl")
include("Jacobian.jl")
include("CodeGeneration.jl")

# FastDifferentiationVisualizationExt overloads them
function make_dot_file end
function draw_dot end
function write_dot end


function test()
    x, graph, _, _ = simple_dominator_graph()

    factor!(graph)
    fedge = edges(graph, 1, 4)[1]
    tmp0 = make_function([value(fedge)], [x])
    dfsimp(x) = tmp0([x])[1]
    x, graph, _, _ = simple_dominator_graph() #x is a new variable so have to make a new Node(x)

    tmp00 = make_function([root(graph, 1)], [x])
    origfsimp(x) = tmp00([x])[1]
    @assert isapprox(FiniteDifferences.central_fdm(5, 1)(origfsimp, 3), dfsimp(3)[1])
end

include("FDTests.jl")


end # module

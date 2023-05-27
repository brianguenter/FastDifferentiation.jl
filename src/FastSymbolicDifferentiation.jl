module FastSymbolicDifferentiation

using TermInterface
using Symbolics: Num, @variables, NoDeriv
import SymbolicUtils
using SymbolicUtils: arguments
using StaticArrays
using SpecialFunctions
using NaNMath
using Statistics
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
include("GraphProcessing.jl")
include("FactorableSubgraph.jl")
include("Factoring.jl")
include("GraphVisualization.jl")

include("FSDTests.jl")
include("Scratchpad.jl")

end # module

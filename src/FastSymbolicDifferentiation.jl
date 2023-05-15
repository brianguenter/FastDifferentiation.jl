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
using Logging
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

"""
Created this because it is not safe to use @infiltrate. If you use `@infiltrate` then you also must have a `using` or `import Infiltrator` 
statement. Unfortunately this triggers the package manager to 
add `Infiltrator` as a direct dependency in the Project.toml file for `FastSymbolicDifferentiation`.

Unfortunately package manager won't remove `Infiltrator` 
as a direct dependency even if all the `@infiltrate` and `using Infiltrate` statements are removed. Because of this `Infiltrator` will 
ship as a dependendency in the released code.
Want `Infiltrator` to only be a dependency for the global Julia
environment, not the FastSymbolicDifferentiation environment.
"""
macro safe_infiltrate()
    :(
        if isdefined(Main, :Infiltrator)
            Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
        end
    )
end

RuntimeGeneratedFunctions.init(@__MODULE__)

const DefaultNodeIndexType = Int64

include("Utilities.jl")
include("BitVectorFunctions.jl")
include("ExpressionGraph.jl")
include("UnspecifiedFunction.jl")
include("PathEdge.jl")
include("DerivativeGraph.jl")
include("GraphProcessing.jl")
include("FactorableSubgraph.jl")
include("Factoring.jl")
include("GraphVisualization.jl")

include("FSDTests.jl")
include("Scratchpad.jl")

# logger = SimpleLogger(open("log1.txt", "w+"))

# global_logger(logger)
# @info "Testing logger"

end # module

module FastSymbolicDifferentiation

using TermInterface
using Symbolics: Num, @variables, NoDeriv
import SymbolicUtils
using SymbolicUtils: arguments
using StaticArrays
using SpecialFunctions
using NaNMath
using TimerOutputs #WARNING:TODO remember to remove this when done benchmarking
using Statistics
using RuntimeGeneratedFunctions


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
export safe_infiltrate

RuntimeGeneratedFunctions.init(@__MODULE__)

const DefaultNodeIndexType = Int64

include("Utilities.jl")
include("BitVectorFunctions.jl")
include("AutomaticDifferentiation.jl")
include("PathEdge.jl")
include("RnToRmGraph.jl")
include("GraphProcessing.jl")
include("FactorableSubgraph.jl")
include("Factoring.jl")
include("Analysis.jl")
include("SphericalHarmonics.jl")
include("GraphVisualization.jl")

include("Tests.jl")
include("Scratchpad.jl")



end # module

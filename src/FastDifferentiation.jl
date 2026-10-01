module FastDifferentiation

# using TermInterface
using StaticArrays
using SpecialFunctions
using NaNMath
# using Statistics
using RuntimeGeneratedFunctions
import Base: iterate
using UUIDs
using SparseArrays
using DataStructures
import DiffRules

module AutomaticDifferentiation
struct NoDeriv
end
export NoDeriv
end #module

RuntimeGeneratedFunctions.init(@__MODULE__)

const DefaultNodeIndexType = Int64

include("Methods.jl") #functions and macros to generate Node specialized methods for all the common arithmetic, trigonometric, etc., operations.
include("Utilities.jl")
include("BitVectorFunctions.jl")
include("ExpressionGraph.jl") #definition of Node type from which FD expression graphs are created
include("DifferentiationRules.jl")
include("PathEdge.jl")  #functions to create and manipulate edges in derivative graphs
include("Conditionals.jl")
include("DerivativeGraph.jl") #functions to compute derivative graph from an expression graph of Node
include("Reverse.jl") #symbolic implementation of conventional reverse automatic differentiation
include("GraphProcessing.jl")
include("FactorableSubgraph.jl")
include("Factoring.jl")
include("Jacobian.jl") #functions to compute jacobians, gradients, hessians, etc.
include("CodeGeneration.jl") #functions to convert expression graphs of Node to executable functions

# FastDifferentiationVisualizationExt overloads them
function make_dot_file end
function draw_dot end
function write_dot end

end # module

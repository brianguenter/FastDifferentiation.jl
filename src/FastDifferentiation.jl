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
module AutomaticDifferentiation
struct NoDeriv
end
export NoDeriv
end #module

const INVARIANTS = true

"""
    @invariant ex msgs...

This macro is used to create invariant test code that is dependent on the global constant `INVARIANTS`. If `INVARIANTS` is false then the test code will not be inserted into the program and there will be no run time overhead. If `INVARIANTS` is true then the code will be inserted. Code that tests invariants tends to increase run time substantially so only set `INVARIANTS` true when you are debugging or testing."""
macro invariant(ex, msgs...)
    if INVARIANTS
        return :(@assert $(esc(ex)) $(esc(msgs)))
    end
end



RuntimeGeneratedFunctions.init(@__MODULE__)

const DefaultNodeIndexType = Int64

include("Methods.jl")
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

include("FDTests.jl")

function test()
    p = make_variables(:p, 21)

    println("NO array zero statement")
    show(make_Expr(p, p, true, true))
    show(make_Expr(p, p, true, false))
    show(make_Expr(p, p, false, true))
    show(make_Expr(p, p, false, false))


    p[21] = 0

    println("shouldn't have an array zero statement but it should have a p[21]= 0 statement")
    show(make_Expr(p, p, true, true))
    println("this should not have an array zero statement nor should have a p[21] = 0 statement")
    show(make_Expr(p, p, true, false))
    println("should not have an array zero statement but should have a p[21] = 0 statement")
    show(make_Expr(p, p, false, true))
    show(make_Expr(p, p, false, false))

    p[20] = 0
    println("this should have an array zero statement should not have p[20]=0 or p[21]=0 statementt")
    show(make_Expr(p, p, true, true))
    println("this should not have an array zero statement should not have p[20]=0 or p[21]=0 statement")
    show(make_Expr(p, p, true, false))
    println("these should both have an array zero creation but should not have p[20]=0 or p[21]=0 statement")
    show(make_Expr(p, p, false, true))
    show(make_Expr(p, p, false, false))
end
export test

end # module

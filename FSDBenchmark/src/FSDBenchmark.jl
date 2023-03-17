module FSDBenchmark
using Symbolics
using FileIO
using BenchmarkTools
using Statistics
using DataFrames
using CSV
using CurveFit
using Plots
using Memoize
using FastSymbolicDifferentiation
using StaticArrays

include("Chebyshev.jl")
include("SphericalHarmonics.jl")
include("StructureFromMotion.jl")
include("LagrangianDynamics.jl")

@variables x, y, z

preprocess_trial(t::BenchmarkTools.Trial, SHOrder::AbstractString) =
    (SHOrder=SHOrder,
        minimum=minimum(t.times),
        median=median(t.times),
        maximum=maximum(t.times),
        allocations=t.allocs,
        memory_estimate=t.memory)
export preprocess_trial

function symbolics_spherical_harmonics(filename, min_order, max_order)
    output = DataFrame()

    for n in min_order:2:max_order
        trial = @benchmark SHDerivatives($n, $x, $y, $z)
        push!(output, preprocess_trial(trial, "$n"))
    end
    CSV.write(filename, output)
    return output
end
export symbolics_spherical_harmonics

function FSD_spherical_harmonics(filename, min_order, max_order)
    output = DataFrame()

    for n in min_order:2:max_order
        trial = @benchmark to_graph($n)
        push!(output, preprocess_trial(trial, "$n"))
    end
    CSV.write(filename, output)
    return output
end
export FSD_spherical_harmonics

function benchmark_spherical_harmonics(min_order, max_order)
    symbolics_spherical_harmonics("Data/SymbolicsSH.csv", min_order, max_order)
    FSD_spherical_harmonics("Data/FSDSH.csv", min_order, max_order)
end
export benchmark_spherical_harmonics


function plot_data(bench1, bench2)
    data1 = CSV.read(bench1, DataFrame)
    data2 = CSV.read(bench2, DataFrame)

    # plot(data1[:, :SHOrder], data1[:, :minimum] / 1e6, ylabel="ms", xlabel="Spherical Harmonic Order")
    # plot!(data2[:, :SHOrder], data2[:, :minimum] / 1e6, ylabel="ms", xlabel="Spherical Harmonic Order")
    plot(data1[:, :SHOrder], data2[:, :minimum] ./ data1[:, :minimum], xlabel="Spherical Harmonic Order", ylabel="Speed of FSD relative to Symbolics",)
end
export plot_data

end # module Benchmarks

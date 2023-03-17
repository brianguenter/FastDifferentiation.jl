module DifferentiationBenchmarks
using Symbolics
using FileIO
using BenchmarkTools
using Statistics
using DataFrames
using CSV
using CurveFit
using Plots

@variables x, y, z

preprocess_trial(t::BenchmarkTools.Trial, SHOrder::AbstractString) =
    (SHOrder=SHOrder,
        minimum=minimum(t.times),
        median=median(t.times),
        maximum=maximum(t.times),
        allocations=t.allocs,
        memory_estimate=t.memory)
export preprocess_trial

function Symbolics_Spherical_Harmonics(filename, min_order, max_order)
    output = DataFrame()

    for n in min_order:2:max_order
        trial = @benchmark SHDerivatives($n, $x, $y, $z)
        push!(output, preprocess_trial(trial, "$n"))
    end
    CSV.write(filename, output)
    return output
end
export Symbolics_Spherical_Harmonics

function FSD_Spherical_Harmonics(filename, min_order, max_order)
    output = DataFrame()

    for n in min_order:2:max_order
        trial = @benchmark to_graph($n)
        push!(output, preprocess_trial(trial, "$n"))
    end
    CSV.write(filename, output)
    return output
end
export FSD_Spherical_Harmonics

function Run_Spherical_Harmonics(min_order, max_order)
    Symbolics_Spherical_Harmonics("Data/SymbolicsSH.csv", min_order, max_order)
    FSD_Spherical_Harmonics("Data/FSDSH.csv", min_order, max_order)
end
export Run_Spherical_Harmonics


function plot_data(bench1, bench2)
    data1 = CSV.read(bench1, DataFrame)
    data2 = CSV.read(bench2, DataFrame)

    # plot(data1[:, :SHOrder], data1[:, :minimum] / 1e6, ylabel="ms", xlabel="Spherical Harmonic Order")
    # plot!(data2[:, :SHOrder], data2[:, :minimum] / 1e6, ylabel="ms", xlabel="Spherical Harmonic Order")
    plot(data1[:, :SHOrder], data1[:, :minimum] ./ data2[:, :minimum], xlabel="Spherical Harmonic Order", ylabel="Relative Speed",)
end
export plot_data

end # module Benchmarks

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

function create_Symbolics_exe(max_l)
    @variables x, y, z

    jac = Symbolics.jacobian(SHFunctions(max_l, x, y, z), [x, y, z]; simplify=true)
    fn1, fn2 = eval.(build_function(jac, [x, y, z]))
    return fn1, fn2
end
export create_Symbolics_exe

function make_FSD_exe(max_l)
    graph, x, y, z = to_graph(max_l)
    return jacobian_function!(graph, Node.([x, y, z]))
end
export make_FSD_exe

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

    for n in min_order:1:max_order
        trial = @benchmark Symbolics.jacobian(fn, [$x, $y, $z], simplify=true) setup = fn = create_Symbolics_SH_functions($n, $x, $y, $z)
        push!(output, preprocess_trial(trial, "$n"))
    end
    CSV.write(filename, output)
    return output
end
export symbolics_spherical_harmonics

function FSD_spherical_harmonics(filename, min_order, max_order)
    output = DataFrame()

    for n in min_order:1:max_order
        trial = @benchmark FastSymbolicDifferentiation.symbolic_jacobian!(gr) setup = gr = to_graph($n)[1] evals = 1
        push!(output, preprocess_trial(trial, "$n"))
    end
    CSV.write(filename, output)
    return output
end
export FSD_spherical_harmonics

function benchmark_spherical_harmonics(min_order, max_order)
    FSD_spherical_harmonics("Data/FSDSH.csv", min_order, max_order)
    symbolics_spherical_harmonics("Data/SymbolicsSH.csv", min_order, max_order)
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

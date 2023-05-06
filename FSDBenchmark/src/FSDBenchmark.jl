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
import FastSymbolicDifferentiation
using FastSymbolicDifferentiation: derivative, jacobian_function!, symbolic_jacobian!, Node, UnspecifiedFunction, codomain_dimension, domain_dimension, function_of, number_of_operations, DerivativeGraph
using StaticArrays
using LaTeXStrings
import LinearAlgebra


const FSD = FastSymbolicDifferentiation

include("Chebyshev.jl")
include("SphericalHarmonics.jl")
include("Transformations.jl")
include("LagrangianDynamics.jl")
include("SimpsonHermite.jl")


@variables x, y, z

const SYMBOLIC = "symbolic"
const EXE = "exe"
const MAKE_FUNCTION = "make_function"

filename(function_name, benchmark_name, min_order, max_order, simplify) = "Data/$(function_name)-$(benchmark_name)-$min_order-$max_order-simplification-$simplify.csv"


function plot_data(bench1, bench2, graph_title, simplify)
    data1 = CSV.read(bench1, DataFrame)
    data2 = CSV.read(bench2, DataFrame)

    graph_title = "Ratio of times, Symbolics/FSD: $graph_title"
    # plot(data1[:, :SHOrder], data1[:, :minimum] / 1e6, ylabel="ms", xlabel="Spherical Harmonic Order")
    # plot!(data2[:, :SHOrder], data2[:, :minimum] / 1e6, ylabel="ms", xlabel="Spherical Harmonic Order")
    p = plot(data1[:, :SHOrder], data2[:, :minimum] ./ data1[:, :minimum], xlabel="Spherical Harmonic Order", ylabel="Ratio", title=graph_title, titlefontsizes=10, legend=false)

    return p
end
export plot_data


function plot_SH_times(min_order, max_order, simplify, graph_title, filename_function, algorithm_name)
    msec_scale = 1e-6
    processes = ["symbolic", "make_function", "exe"]
    for process in processes
        data = CSV.read(filename_function(SH_NAME, process, min_order, max_order, simplify), DataFrame)
        plot!(data[:, :minimum] * msec_scale, label=process * algorithm_name)
    end

    plot!(xlabel="Spherical Harmonic Order", yscale=:log10, ylabel="log₁₀(time),msecs", title=graph_title, legend_position=:topleft)
end
export plot_SH_times

function plot_SH_symbolic_time(min_order, max_order, simplify)
    plot_data(
        FSD_filename(SH_NAME, "symbolic", min_order, max_order, simplify),
        Symbolics_filename(SH_NAME, "symbolic", min_order, max_order, simplify),
        "symbolic processing",
        simplify)
    plot!(xticks=min_order:2:max_order)
end
export plot_SH_symbolic_time

function plot_SH_exe_time(min_order, max_order, simplify)
    plot_data(
        FSD_filename(SH_NAME, EXE, min_order, max_order, simplify),
        Symbolics_filename(SH_NAME, EXE, min_order, max_order, simplify),
        "execution",
        simplify)
    plot!(xticks=min_order:2:max_order)
end
export plot_SH_exe_time

"""plot that shows how FSD jacobian is close to optimal for SH because number of operations of Jacobian is a fixed constant (roughly 2.5) times the number of operations in the original function"""
function plot_SH_FSD_graph_vs_jacobian_size(min_order, max_order)
    funcs = [to_graph(x)[1] for x in min_order:max_order]
    derivs = [symbolic_jacobian!(x) for x in funcs]
    ratio = number_of_operations.(derivs) ./ number_of_operations.(roots.(funcs))
    plot(min_order:max_order, ratio, ylabel="Ratio of operations", title=L"\frac{operations(jacobian(f))}{operations(f)}", titlefontsizes=10, xlabel="Spherical Harmonic order", legend=false)
end
export plot_SH_FSD_graph_vs_jacobian_size

function plot_SH_make_function_time(min_order, max_order, simplify)
    plot_data(
        FSD_filename(SH_NAME, MAKE_FUNCTION, min_order, max_order, simplify),
        Symbolics_filename(SH_NAME, MAKE_FUNCTION, min_order, max_order, simplify),
        "make_function",
        simplify)
    plot!(xticks=min_order:2:max_order)
end
export plot_SH_make_function_time

extract_info(model_size, benchmark_timing) = [model_size, minimum(benchmark_timing).time, median(benchmark_timing).time, maximum(benchmark_timing).time, benchmark_timing.allocs, benchmark_timing.memory]

function write_data(data, model_function, benchmark_function, min_model_size, max_model_size, simplify)
    CSV.write(
        filename(nameof(model_function), nameof(benchmark_function), min_model_size, max_model_size, simplify),
        data)
end

#For FSD functions model_function returns a DerivativeGraph. For Symbolics functions model_function returns a function (may be an array in the function has multiple outputs) as well as a vector of input variables.

fsd_symbolic_benchmark(model_function::Function, model_size; simplify=false) = @benchmark symbolic_jacobian!(gr) setup = gr = $model_function($model_size) evals = 1
fsd_make_benchmark(model_function::Function, model_size; simplify=false) = @benchmark jacobian_function!(graph) setup = (graph = $model_function($model_size)) evals = 1
function fsd_exe_benchmark(model_function::Function, model_size, simplify=false)
    graph = model_function(model_size)
    exe = jacobian_function!(graph, FSD.variables(graph), in_place=true)

    input = rand(domain_dimension(graph))
    output = rand(codomain_dimension(graph), domain_dimension(graph))
    return @benchmark $exe($input..., $output)
end


function symbolics_symbolic_benchmark(model_function, model_size; simplify=false)
    tmp = model_function(model_size)
    println(Symbolics.jacobian(tmp[1], tmp[2]; simplify=simplify))
    @benchmark Symbolics.jacobian(tmp[1], tmp[2]; simplify=$simplify) setup = (tmp = $model_function($model_size))
end

function symbolics_exe_benchmark(model_function, model_size; simplify=false)
    model, vars = model_function(model_size)
    jac = Symbolics.jacobian(model, vars; simplify=simplify)
    out_of_place, in_place = build_function(jac, vars; expression=Val{false})
    tmp_matrix = out_of_place(rand(length(vars)))

    return @benchmark $in_place($tmp_matrix, rand(length($vars)))
end

function symbolics_make_benchmark(model_function, model_size; simplify=false)
    model, vars = model_function(model_size)
    jac = Symbolics.jacobian(model, vars; simplify=simplify)
    return @benchmark build_function($jac, $vars; expression=Val{false})
end

function single_benchmark(model_function::Function, benchmark_function::Function, model_range, simplify=false)
    # @benchmark Symbolics.jacobian(fn, [$x, $y, $z], simplify=$simplify) setup = gr = evals = 1
    data = DataFrame(model_size=Int64[], minimum=Float64[], median=Float64[], maximum=Float64[], allocations=Int64[], memory_estimate=Int64[])

    for model_size in model_range
        time = benchmark_function(model_function, model_size)
        push!(data, extract_info(model_size, time))
    end
    write_data.([data], model_function, benchmark_function, minimum(model_range), maximum(model_range), simplify)
end
export single_benchmark

function benchmark(models, sizes, benchmarks; simplify::Bool=false)
    for (model, size_range) in zip(models, sizes)
        for bench in benchmarks
            single_benchmark(model, bench, size_range, simplify)
        end
    end
end
export benchmark

benchmark_sizes = [5:1:25, 5:5:40]
FSD_params() = ([FSD_spherical_harmonics, FSD_chebyshev], benchmark_sizes, [fsd_symbolic_benchmark, fsd_make_benchmark, fsd_exe_benchmark])
export FSD_params
Symbolics_params() = ([Symbolics_spherical_harmonics, Symbolics_chebyshev], benchmark_sizes, [symbolics_symbolic_benchmark, symbolics_make_benchmark, symbolics_exe_benchmark])
export Symbolics_params

end # module Benchmarks

module FSDBenchmark
using Symbolics
using FileIO
using BenchmarkTools
using Statistics

using CurveFit
using Plots
using Memoize
import FastSymbolicDifferentiation
using StaticArrays
using FastSymbolicDifferentiation: derivative, jacobian_function!, symbolic_jacobian!, Node, UnspecifiedFunction, codomain_dimension, domain_dimension, function_of, number_of_operations, DerivativeGraph
using LaTeXStrings
import LinearAlgebra

using DataFrames
using CSV

filename(model_function, package, benchmark, min_model_size, max_model_size, simplify) = "Data/" * join([nameof(model_function), nameof(typeof(package)), nameof(typeof(benchmark))], "_") * "_$(min_model_size)_$(max_model_size)_simplify_$simplify.csv"
export filename

function write_data(data::DataFrame, model_function, package, benchmark, min_model_size, max_model_size, simplify)
    CSV.write(
        filename(model_function, package, benchmark, min_model_size, max_model_size, simplify),
        data)
end
export write_data

const FSD = FastSymbolicDifferentiation

include("Types.jl")
include("Chebyshev.jl")
include("SphericalHarmonics.jl")
include("Transformations.jl")
include("LagrangianDynamics.jl")
include("SimpsonHermite.jl")


@variables x, y, z



const SYMBOLIC = "symbolic"
const EXE = "exe"
const MAKE_FUNCTION = "make_function"



extract_info(model_size, benchmark_timing) = Any[Int64(model_size), Float64(minimum(benchmark_timing).time), Float64(median(benchmark_timing).time), Float64(maximum(benchmark_timing).time), Int64(benchmark_timing.allocs), Int64(benchmark_timing.memory)]
export extract_info

"""plot that shows how FSD jacobian is close to optimal for SH because number of operations of Jacobian is a fixed constant (roughly 2.5) times the number of operations in the original function"""
function plot_SH_FSD_graph_vs_jacobian_size(min_order, max_order)
    funcs = [to_graph(x)[1] for x in min_order:max_order]
    derivs = [symbolic_jacobian!(x) for x in funcs]
    ratio = number_of_operations.(derivs) ./ number_of_operations.(roots.(funcs))
    plot(min_order:max_order, ratio, ylabel="Ratio of operations", title=L"\frac{operations(jacobian(f))}{operations(f)}", titlefontsizes=10, xlabel="Spherical Harmonic order", legend=false)
end
export plot_SH_FSD_graph_vs_jacobian_size



#Benchmark code for FSD

run_benchmark(model_function::Function, model_size, package::FastSymbolic, ::Symbolic; simplify=false) = @benchmark symbolic_jacobian!(gr) setup = gr = $model_function($package, $model_size) evals = 1

run_benchmark(model_function::Function, model_size, package::FastSymbolic, ::MakeFunction; simplify=false) = @benchmark jacobian_function!(graph) setup = (graph = $model_function($package, $model_size)) evals = 1

function run_benchmark(model_function::Function, model_size, package::FastSymbolic, ::Exe, simplify=false)
    graph = model_function(package, model_size)
    exe = jacobian_function!(graph, FSD.variables(graph), in_place=true)

    input = rand(domain_dimension(graph))
    output = rand(codomain_dimension(graph), domain_dimension(graph))
    return @benchmark $exe($input..., $output)
end

#Benchmark code for Symbolics

function run_benchmark(model_function, model_size, package::JuliaSymbolics, ::Symbolic; simplify=false)
    tmp = model_function(package, model_size)
    @benchmark Symbolics.jacobian(tmp[1], tmp[2]; simplify=$simplify) setup = (tmp = $model_function($package, $model_size))
end

function run_benchmark(model_function, model_size, package::JuliaSymbolics, ::Exe; simplify=false)
    model, vars = model_function(package, model_size)
    jac = Symbolics.jacobian(model, vars; simplify=simplify)
    out_of_place, in_place = build_function(jac, vars; expression=Val{false})
    tmp_matrix = out_of_place(rand(length(vars))) #generate a matrix of the correct size

    return @benchmark $in_place($tmp_matrix, rand(length($vars)))
end

function run_benchmark(model_function, model_size, package::JuliaSymbolics, ::MakeFunction; simplify=false)
    model, vars = model_function(package, model_size)
    jac = Symbolics.jacobian(model, vars; simplify=simplify)
    return @benchmark build_function($jac, $vars; expression=Val{false})
end
export run_benchmark

make_data() = DataFrame(model_size=Int64[], minimum=Float64[], median=Float64[], maximum=Float64[], allocations=Int64[], memory_estimate=Int64[])

function single_benchmark(model_function::Function, model_range, package::AbstractPackage, benchmark::AbstractBenchmark, simplify=false)
    data = make_data()

    for model_size in model_range
        timing = run_benchmark(model_function, model_size, package, benchmark)
        push!(data, extract_info(model_size, timing))
    end

    write_data(data, model_function, package, benchmark, minimum(model_range), maximum(model_range), simplify)
end
export single_benchmark

# function benchmark(models, sizes, package::AbstractPackage, benchmarks::AbstractVector{AbstractBenchmark}; simplify::Bool=false)
#     for (model, size_range) in zip(models, sizes)
#         for bench in benchmarks
#             single_benchmark(model, size_range, package, bench, simplify)
#         end
#     end
# end

benchmark_sizes() = [5:1:25, 5:5:50]
model_functions() = [spherical_harmonics, chebyshev]
benchmark_types() = [Symbolic(), Exe(), MakeFunction()]
export benchmark_types
params() = (model_functions(), benchmark_sizes())
export params

benchmark_package(package, range, model_function; simplify=false) = single_benchmark.(Ref(model_function), Ref(range), Ref(package), benchmark_types(), simplify)


benchmark_package(package; simplify=false) = benchmark_package.(Ref(package), benchmark_sizes(), model_functions(), simplify=simplify)


function plot_data(model_function, bench1, graph_title::AbstractString, xlabel::AbstractString, simplify)
    fname1 = filename(model_function, FastSymbolic(), bench1, extrema(benchmark_sizes()[findfirst(x -> x == model_function, model_functions())])..., simplify)
    fname2 = filename(model_function, JuliaSymbolics(), bench1, extrema(benchmark_sizes()[findfirst(x -> x == model_function, model_functions())])..., simplify)
    data1 = CSV.read(fname1, DataFrame)
    data2 = CSV.read(fname2, DataFrame)

    graph_title = "Ratio of times, Symbolics/FSD: $graph_title"
    println("here")
    # plot(data1[:, :SHOrder], data1[:, :minimum] / 1e6, ylabel="ms", xlabel="Spherical Harmonic Order")
    # plot!(data2[:, :SHOrder], data2[:, :minimum] / 1e6, ylabel="ms", xlabel="Spherical Harmonic Order")
    p = plot(data1[:, :model_size], data2[:, :minimum] ./ data1[:, :minimum], xlabel=xlabel, ylabel="Ratio", title=graph_title, titlefontsizes=10, legend=false, marker=:circle)

    return p
end
export plot_data

function test_Symbolics_limit()
    for i in 20:25
        try
            model, vars = spherical_harmonics(JuliaSymbolics(), i)
            jac = Symbolics.jacobian(model, vars; simplify=false)
            out_of_place, in_place = build_function(jac, vars; expression=Val{false})
            println("finished size $i")
        catch exc
            println("failed at size $i $exc")
            break
        end
    end
end
export test_Symbolics_limit


end # module Benchmarks

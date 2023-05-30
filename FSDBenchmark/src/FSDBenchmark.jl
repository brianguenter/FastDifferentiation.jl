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
using FastSymbolicDifferentiation: derivative, jacobian_function!, symbolic_jacobian!, Node, codomain_dimension, domain_dimension, number_of_operations, DerivativeGraph
using LaTeXStrings
import LinearAlgebra

using DataFrames
using CSV

filename(model_function, package, benchmark, min_model_size, max_model_size, simplify) = DATA_DIR * join([nameof(model_function), nameof(typeof(package)), nameof(typeof(benchmark))], "_") * "_$(min_model_size)_$(max_model_size)_simplify_$simplify.csv"
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
include("SimpsonHermite.jl")


@variables x, y, z


const DATA_DIR = "Data/"
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

## Benchmark code for FSD
run_benchmark(model_function::Function, model_size, package::FastSymbolic, ::Symbolic; simplify=false) = @benchmark symbolic_jacobian!(gr) setup = gr = $model_function($package, $model_size) evals = 1

run_benchmark(model_function::Function, model_size, package::FastSymbolic, ::MakeFunction; simplify=false) = @benchmark jacobian_function!(graph) setup = (graph = $model_function($package, $model_size)) evals = 1

function run_benchmark(model_function::Function, model_size, package::FastSymbolic, ::Exe, simplify=false)
    graph = model_function(package, model_size)
    exe = jacobian_function!(graph, FSD.variables(graph), in_place=true)

    input = rand(domain_dimension(graph))
    output = rand(codomain_dimension(graph), domain_dimension(graph))
    return @benchmark $exe($input..., $output)
end
## end benchmark functions for FSD


## Benchmark code for Symbolics
function run_benchmark(model_function, model_size, package::JuliaSymbolics, ::Symbolic; simplify=false)
    tmp = model_function(package, model_size)
    @benchmark Symbolics.jacobian(tmp[1], tmp[2]; simplify=$simplify) setup = (tmp = $model_function($package, $model_size))
end

function run_benchmark(model_function, model_size, package::JuliaSymbolics, ::Exe; simplify=false)
    @info "creating executable"
    model, vars = model_function(package, model_size)
    jac = Symbolics.jacobian(model, vars; simplify=simplify)
    out_of_place, in_place = build_function(jac, vars; expression=Val{false})
    tmp_matrix = Matrix{Float64}(undef, length(model), length(vars)) #generate a matrix of the correct size.
    @info "done creating executable, starting Exe benchmark"

    return @benchmark $in_place($tmp_matrix, rand(length($vars)))
end

function run_benchmark(model_function, model_size, package::JuliaSymbolics, ::MakeFunction; simplify=false)
    model, vars = model_function(package, model_size)
    jac = Symbolics.jacobian(model, vars; simplify=simplify)
    return @benchmark build_function($jac, $vars; expression=Val{false})
end
export run_benchmark
## End benchmark functions for Symbolics.jl

## Generic benchmark functions
make_data() = DataFrame(model_size=Int64[], minimum=Float64[], median=Float64[], maximum=Float64[], allocations=Int64[], memory_estimate=Int64[])

function single_benchmark(model_function::Function, model_range, package::AbstractPackage, benchmark::AbstractBenchmark, simplify=false)
    data = make_data()
    FSD.clear_cache() #clear cache otherwise memory usage keeps rising monotonically. Wouldn't expect this to be a problem until functions become gigantic, billions of nodes.
    for model_size in model_range
        @info "Starting model size $model_size"
        timing = run_benchmark(model_function, model_size, package, benchmark, simplify=simplify)
        push!(data, extract_info(model_size, timing))
        @info "Finished for model size $model_size"

        #incrementally write benchmarks out in case something crashes or benchmarks take too long to complete the entire run.
        CSV.write(
            filename(model_function, package, benchmark, minimum(model_range), model_size, simplify),
            data)
        if model_size != minimum(model_range)
            rm(filename(model_function, package, benchmark, minimum(model_range), model_size - 1, simplify)) #delete the previous file to avoid cluttering the directory
        end
    end
end
export single_benchmark

benchmark_sizes() = [5:1:25, 5:1:30]
export benchmark_sizes
model_functions() = [spherical_harmonics, chebyshev]
benchmark_types() = [Symbolic(), Exe(), MakeFunction()]
export benchmark_types
params() = (model_functions(), benchmark_sizes())
export params

"""More specialized version of `benchmark_package` that allows you to choose which range and model function to use"""
benchmark_package(package, range, model_function; simplify=false) = single_benchmark.(Ref(model_function), Ref(range), Ref(package), benchmark_types(), simplify)

"""For the package defined by `package` run all benchmarks defined by `benchmark_types()` on the functions defined by `model_functions()` for the benchmark sizes defined by `benchmark_sizes()`. The arrays returned by `benchmark_sizes()` and `model_sizes()` are aligned. For example the `spherical_harmonics` model function which is model_functions()[1] will be run with model sizes `benchmark_sizes()[1]`.

The legal package and benchmark types are defined in Types.jl. `package` can be one of `Julia_Symbolics(), FastSymbolic()`, where `JuliaSymbolics` runs the benchmarks using `Symbolics.jl`.

If `simplify=false` then simplification will not be used in the `Symbolics.jl` benchmarks. If it is true then simplified will be used. The latter can significantly increase run time of the benchmarks.

Example:
```
benchmark_package(JuliaSymbolics())
```
will run the benchmarks`[Symbolic(), Exe(), MakeFunction()]` on the model functions `[spherical_harmonics, chebyshev]`. The `spherical_harmonics` model function will be called with model sizes 5:1:25 and the `chebyshev` model function will be called with model sizes 5:1:30.
"""
benchmark_package(package; simplify=false) = benchmark_package.(Ref(package), benchmark_sizes(), model_functions(), simplify=simplify)
export benchmark_package

"""Run all benchmarks for all packages"""
benchmark_all(simplify::Bool) = benchmark_package.([FastSymbolic(), JuliaSymbolics()]; simplify=simplify)
export benchmark_all

function plot_data(model_function, bench1, simplify)
    benchmark_name = nameof(typeof(bench1))
    fname1 = filename(model_function, FastSymbolic(), bench1, extrema(benchmark_sizes()[findfirst(x -> x == model_function, model_functions())])..., false)
    fname2 = filename(model_function, JuliaSymbolics(), bench1, extrema(benchmark_sizes()[findfirst(x -> x == model_function, model_functions())])..., simplify)
    data1 = CSV.read(fname1, DataFrame)
    data2 = CSV.read(fname2, DataFrame)
    println(data2)
    graph_title = "Time ratio, Symbolics/FSD \n$(nameof(model_function)) \n$benchmark_name benchmark, simplify = $simplify"

    #find first missing value
    last_good = findfirst(x -> x === missing, data2[:, :minimum])
    if last_good === nothing
        last_good = size(data2)[1]
    end

    println("lastgood $last_good")
    ratio = [data2[i, :minimum] === missing ? missing : data2[i, :minimum] / data1[i, :minimum] for i in 1:size(data2)[1]]
    println("ratio $ratio")
    p = plot(data1[1:last_good, :model_size], ratio, xlabel="$(nameof(model_function)) model size", ylabel="Ratio", title=graph_title, titlefontsizes=10, legend=false, marker=:circle, yaxis=:log, minorticks=5)

    return p
end
export plot_data

function publication_benchmarks(simplify::Bool, model_functions::AbstractVector=model_functions(), benchmarks::AbstractVector=benchmark_types(), run_benchmarks::Bool=true)
    if run_benchmarks
        benchmark_all(simplify)
    end

    for bench in benchmarks
        for model in model_functions
            bench_type = typeof(bench)
            savefig(plot_data(model, bench, simplify),
                "$(DATA_DIR)figure_$(model)_$(bench_type)_simplify_$(simplify).svg"
            )
        end
    end
end
export publication_benchmarks

end # module Benchmarks

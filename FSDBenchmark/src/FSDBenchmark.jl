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
include("Transformations.jl")
include("LagrangianDynamics.jl")
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
    @info "creating executable"
    model, vars = model_function(package, model_size)
    jac = Symbolics.jacobian(model, vars; simplify=simplify)
    out_of_place, in_place = build_function(jac, vars; expression=Val{false})
    tmp_matrix = Matrix{Float64}(undef, length(model), length(vars)) #generate a matrix of the correct size.
    @info "done creating executable, starting Exe benchmark"

    return @benchmark $in_place($tmp_matrix, rand(length($vars)))
end

function run_benchmark(model_function, model_size, package::JuliaSymbolics, ::MakeFunction; simplify=false)
    @info "computing symbolic jacobian"
    model, vars = model_function(package, model_size)
    jac = Symbolics.jacobian(model, vars; simplify=simplify)
    @info "done computing jacobian beginning MakeFunction benchmark"
    return @benchmark out_of_place, in_place = build_function($jac, $vars; expression=Val{false})
end
export run_benchmark

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
benchmark_types() = [Symbolic(), MakeFunction(), Exe()]
export benchmark_types
params() = (model_functions(), benchmark_sizes())
export params

benchmark_package(package, range, model_function; simplify=false) = single_benchmark.(Ref(model_function), Ref(range), Ref(package), benchmark_types(), simplify)


benchmark_package(package; simplify=false) = benchmark_package.(Ref(package), benchmark_sizes(), model_functions(), simplify=simplify)
export benchmark_package

benchmark_all(simplify::Bool) = benchmark_package.([FastSymbolic(), JuliaSymbolics()]; simplify=simplify)
export benchmark_all

function plot_data(model_function, bench1, graph_title::AbstractString, xlabel::AbstractString, simplify)
    fname1 = filename(model_function, FastSymbolic(), bench1, extrema(benchmark_sizes()[findfirst(x -> x == model_function, model_functions())])..., simplify)
    fname2 = filename(model_function, JuliaSymbolics(), bench1, extrema(benchmark_sizes()[findfirst(x -> x == model_function, model_functions())])..., simplify)
    data1 = CSV.read(fname1, DataFrame)
    data2 = CSV.read(fname2, DataFrame)
    println(data2)
    graph_title = "Time ratio, Symbolics/FSD: $graph_title"

    #find first missing value
    last_good = findfirst(x -> x === missing, data2[:, :minimum])
    if last_good === nothing
        last_good = size(data2)[1]
    end

    println("lastgood $last_good")
    ratio = [data2[i, :minimum] === missing ? missing : data2[i, :minimum] / data1[i, :minimum] for i in 1:size(data2)[1]]
    println("ratio $ratio")
    p = plot(data1[1:last_good, :model_size], ratio, xlabel=xlabel, ylabel="Ratio", title=graph_title, titlefontsizes=10, legend=false, marker=:circle)

    return p
end
export plot_data

function publication_benchmarks(simplify::Bool, run_benchmarks=true)
    if run_benchmarks
        benchmark_all(simplify)
    end

    for bench in benchmark_types()
        for model in model_functions()
            bench_type = typeof(bench)
            savefig(plot_data(model, bench, "$bench_type\n$model", "$model order", simplify),
                "$(DATA_DIR)figure_$(model)_$(bench_type).svg"
            )
        end
    end
end
export publication_benchmarks

function test_Symbolics_limit()
    model, vars = spherical_harmonics(JuliaSymbolics(), 20)
    jac = Symbolics.jacobian(model, vars; simplify=false)
    out_of_place, in_place = build_function(jac, vars; expression=Val{false})
    tmp_matrix = out_of_place(rand(length(vars))) #generate a matrix of the correct size

    in_place(tmp_matrix, rand(length(vars)))

    # for i in 5:1:50
    #     try
    #         model, vars = chebyshev(JuliaSymbolics(), i)
    #         jac = Symbolics.jacobian(model, vars; simplify=false)
    #         build_function(jac, vars; expression=Val{false})
    #         # out_of_place, in_place = build_function(jac, vars; expression=Val{false})
    #         println("finished size $i")
    #     catch exc
    #         println("failed at size $i $exc")
    #         break
    #     end
    # end
end
export test_Symbolics_limit


end # module Benchmarks

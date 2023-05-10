
abstract type AbstractPackage end
struct JuliaSymbolics <: AbstractPackage end
struct FastSymbolic <: AbstractPackage end
export JuliaSymbolics, FastSymbolic

abstract type AbstractBenchmark end
struct Symbolic <: AbstractBenchmark end
struct Exe <: AbstractBenchmark end
struct MakeFunction <: AbstractBenchmark end
export Symbolic, Exe, MakeFunction
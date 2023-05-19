This subpackage is used to run benchmarks comparing `FastSymbolicDifferentiation.jl` to `Symbolics.jl`. To run the benchmarks:

```
using FSDBenchmark

benchmark_all(false)
```

This will run all the benchmarks on all the model functions for both `FastSymbolicDifferentiation.jl` and `Symbolics.jl`. The argument to `benchmark_all` is a boolean flag turning simplification on or off in Symbolics.jl. If you set it to true the executable results for `Symbolics.jl` will be better because more simplification will be done (in the case of the `chebyshev` model function much better) but the benchmarks may take much longer to run.

There are other benchmark functions in `FastSymbolicDifferentiation.jl` that allow more fine grained control of which functions, packages, etc., will be included in the benchmarks. See the source for details.
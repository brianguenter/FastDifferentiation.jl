## Future work
The top priority work item is reducing memory allocations and improving performance of the preprocessing step.

The **FD** algorithm is fast enough to preprocess expression graphs with ≈10⁵ operations in approximately 1 minute on a modern laptop. This is a one time step to generate compiled derivative functions which, of course, take far less than 1 minute to execute. Preprocessing time scales non-linearly with expression size so smaller graphs take much less time.

However, LLVM compile time can be significant at this scale. For expressions this size and larger [DynamicExpressions.jl](https://github.com/SymbolicML/DynamicExpressions.jl) might be a better tradeoff between compile and execution time. This would be especially useful when your function is changing frequently so compilation overhead cannot be amortized across many derivative evaluations.

**FD** cannot differentiate expressions with conditionals involving **FD** variables. The algorithm can be extended to allow this but symbolic processing can scale exponentially with the nesting depth of conditionals. For small nesting depths this might be acceptable so a future version of FD might support limited nested conditionals. 

However, a better approach might be to use FD as a processing step in a tracing JIT compiler, applying **FD** to the basic blocks detected and compiled by the JIT. These basic blocks do not have branches. Many programs could be differentiated competely automatically by this method. I'm not a compiler expert so it is unlikely I will do this by myself. But contact me if *you* are a compiler expert and want a cool project to work on.
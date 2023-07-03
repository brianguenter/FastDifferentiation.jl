## Future Work
The two biggest limitations of **FD** are no support for conditionals involving **FD** variables and code size that scales with number of operations. 

The first can be dealt with by modifying `make_function`. The new generated code would have three distinct steps:
1. Evaluate conditionals that involve **FD** variables and generate a bit vector of the boolean values.
2. Index into a dictionary with the bit vector to see if this combination of conditionals has been seen before. If it has use the previously generated executable.
3. If necessary generate a new executable where the conditional values are all known. Generate the new expression graph and compute and complile the derivative code and cache it.

The size of the dictionary can be allowed to grow arbitrarily or an LRU type scheme can be used to discard the least recently used executable from the cahce.

If the set of conditials has good temporal locality then it shouldn't be necessary to recompile very often. The amortized cost should be good.

When expression graphs get larger than about 10^5 nodes LLVM compilation time rises rapidly and can be extremely long. Large programs cannot fit in the smallest and fastest caches, so performance will be limited by cache memory bandwidth, because each instruction is executed only once. In the worst case the runtime generated function wouldn't even fit in the L3 cache and performance would be limited by main memory bandwidth.

 A method for dealing with excessive code size is loop rerolling. **FD** essentially unrolls all loops so vector operations become scalar operations. It should be possible to recognize certain common patterns of unrolling and to undo them. Undoing the unrolling all the way back to the original source expression isn't necessary so long as the code size can be substantially reduced.

The most common pattern that causes code expansion is tensor contraction, which occurs in matrix-vector and matrix-matrix operations. These patterns are simple and should be easily recognized and rerolled. The original expression graph can be annotated with metadata that makes rerolling easier. For example:

```julia
julia> A = make_variables(:a,2,2)
2×2 Matrix{FastDifferentiation.Node}:
 a1_1  a1_2
 a2_1  a2_2

julia> b = make_variables(:b,2)
2-element Vector{FastDifferentiation.Node}:
 b1
 b2

julia> jacobian(cos.(A*b),vcat(vec(A),b))
2×6 Matrix{FastDifferentiation.Node}:
 (-(sin(((a1_1 * b1) + (a1_2 * b2)))) * b1)                                         0.0  …  (-(sin(((a1_1 * b1) + (a1_2 * b2)))) * a1_2)
                                        0.0  (-(sin(((a2_1 * b1) + (a2_2 * b2)))) * b1)     (-(sin(((a2_1 * b1) + (a2_2 * b2)))) * a2_2)

julia> A
2×2 Matrix{FastDifferentiation.Node}:
 a1_1  a1_2
 a2_1  a2_2

julia> b
2-element Vector{FastDifferentiation.Node}:
 b1
 b2

julia> A*b
2-element Vector{Any}:
 ((a1_1 * b1) + (a1_2 * b2))
 ((a2_1 * b1) + (a2_2 * b2))

julia> cos.(A*b)
2-element Vector{FastDifferentiation.Node{typeof(cos), 1}}:
 cos(((a1_1 * b1) + (a1_2 * b2)))
 cos(((a2_1 * b1) + (a2_2 * b2)))

julia> jacobian(ans,vcat(vec(A),b))
2×6 Matrix{FastDifferentiation.Node}:
 (-(sin(((a1_1 * b1) + (a1_2 * b2)))) * b1)  0.0                      …                     (-(sin(((a1_1 * b1) + (a1_2 * b2)))) * a1_2)
0.0                                         (-(sin(((a2_1 * b1) + (a2_2 * b2)))) * b1)      (-(sin(((a2_1 * b1) + (a2_2 * b2)))) * a2_2)
```

Because the variable indices are carried in the variable names it should be relatively easy to spot tensor contraction sequences like this `(a1_1 * b1) + (a1_2 * b2))` and replace them with a tensor contraction operator on matrix elements.


Another possibility is to use [DynamicExpressions.jl](https://github.com/SymbolicML/DynamicExpressions.jl) instead of LLVM compilation of large runtime generated programs. When graph size goes above 10^5 nodes LLVM compilation time can increase dramatically, but DynamicExpressions can be quite fast at even larger scales and have reasonable performance. This would be especially useful when the function is changing frequently so compilation overhead cannot be amortized across many derivative evaluations.

Some hybrid of loop rerolling and DynamicExpressions may make it possible to scale **FD** to expressions several orders of magnitude larger than the current practical limit.

In the short term reducing memory allocations during the derivative computation stage should substantially improve performance.



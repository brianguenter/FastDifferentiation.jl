# Limitations
**FD** does not support expressions with conditionals on **FD** variables. For example, you can do this:
```
julia> f(a,b,c) = a< 1.0 ? cos(b) : sin(c)
f (generic function with 2 methods)

julia> f(0.0,x,y)
cos(x)

julia> f(1.0,x,y)
sin(y)
```
but you can't do this:
```
julia> f(a,b) = a < b ? cos(a) : sin(b)
f (generic function with 2 methods)

julia> f(x,y)
ERROR: MethodError: no method matching isless(::FastDifferentiation.Node{Symbol, 0}, ::FastDifferentiation.Node{Symbol, 0})

Closest candidates are:
  isless(::Any, ::DataValues.DataValue{Union{}})
   @ DataValues ~/.julia/packages/DataValues/N7oeL/src/scalar/core.jl:291
  isless(::S, ::DataValues.DataValue{T}) where {S, T}
   @ DataValues ~/.julia/packages/DataValues/N7oeL/src/scalar/core.jl:285
  isless(::DataValues.DataValue{Union{}}, ::Any)
   @ DataValues ~/.julia/packages/DataValues/N7oeL/src/scalar/core.jl:293
  ...
```
This is because the call `f(x,y)` creates an expression graph. At graph creation time the **FD** variables `x,y` are unevaluated variables with no specific value so they cannot be compared with any other value.

The algorithm can be extended to work with conditionals applied to **FD** variables but the processing time and graph size may grow exponentially with conditional nesting depth. A future version may allow for limited conditional nesting. See [Future Work](@ref) for a potential long term solution to this problem.

The preprocessing/compilation step for expressions graphs with more than ≈10⁵ operations may take a minute or more (processing time increase non-linearly with number of expressions so smaller expression have much shorter preprocessing times). This is due to two factors. 

The current code is not memory efficient - it allocates much more than necessary which makes it slower than it should be. Future versions will be more memory efficient.

The code uses BitVector for tracking reachability of function roots and variable nodes. This seemed like a good idea when I began and thought **FD** would only be practical for modest size expressions (<10⁴ operations). But, it scaled better than expected and for larger graphs the memory overhead of the BitVector representation is significant. Using Set instead of BitVector for larger graphs should significantly reduce symbolic processing time.

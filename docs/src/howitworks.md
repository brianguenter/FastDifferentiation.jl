# How it works
The **FD** differentiation algorithm is related to the [D*](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) algorithm but is asymptotically faster so it works on much larger expression graphs. The new algorithms used in **FD** will be described in a soon to be written paper.

**FD** transforms the input expression graph into a derivative graph[^a], and then factors this graph to generate an efficient expression for the derivative. This is fundamentally different from forward and reverse automatic differentiation. 

The efficiency of **FD** comes from analysis of the graph structure of the function rather than sophisticated algebraic simplification rules. By default **FD** applies only these algebraic simplications[^c] to expressions:
* `x×0=>0`
* `x×1=>x`
* `x/1=>x`
* `x+0=>x`
* `c₁×c₂=>c₃` for `c₁,c₂,c₃` constants
* `c₁+c₂=>c₃` for `c₁,c₂,c₃` constants
* `c₁×(c₂×x)` => `(c₁×c₂)×x`  for `c₁,c₂` constants

These rules are generally safe in the sense of obeying IEEE floating point arithmetic rules. However if the runtime value of `x` happens to be NaN or Inf the **FD** expression `x*0` will identically return 0, because it will have been rewritten to 0 by the simplification rules. The expected IEEE result is NaN.
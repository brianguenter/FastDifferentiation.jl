# How it works
**FD** is a domain specific language (DSL) embedded in Julia. **FD** defines a custom `Number` type and nd overloads all the mathematical operators in Base to apply to this new number type. You create **FD** numbers using either [`@variables`](@ref) or [`make_variable`](@ref).

Mathematical operations on **FD** numbers create a graph representing the mathematical expression rather than immediately returning a floating point value. For example, in this code fragment 
```julia
@variables x y
f(a,b)= cos(a)*sin(b)

myexpr = f(x,y)
```
`myexpr` contains a graph representation of the `cos(x)*sin(y)` where `x,y` are **FD** numbers. 

For the most part there is no difference between using **FD** numbers and the base number types, `Float64`, `Int64`, etc. You define your Julia function as you normally world and then call it with **FD** numbers as inputs; the return value will be a graph representing the expression your Julia function computes.  

The **FD** differentiation functions, `jacobian`, `hessian`, etc., take **FD** expression graphs as inputs and return **FD** expression graphs. To turn this into executable Julia code you pass an **FD** expression graph as an argument to `make_function`.

All the **FD** differntiation functions use derivative graph factorization[^2] to compute derivatives. The **FD** differentiation algorithm is related to the [D*](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) algorithm but is asymptotically faster so it works on much larger expression graphs. The new algorithms used in **FD** will be described in a soon to be written paper. **FD** automatic differentiaion is fundamentally different from forward and reverse automatic differentiation. 

The efficiency of **FD** comes from analysis of the graph structure of the function rather than sophisticated algebraic simplification rules. By default **FD** applies only these algebraic simplications[^1] to expressions:
* `x×0=>0`
* `x×1=>x`
* `x/1=>x`
* `x+0=>x`
* `c₁×c₂=>c₃` for `c₁,c₂,c₃` constants
* `c₁+c₂=>c₃` for `c₁,c₂,c₃` constants
* `c₁×(c₂×x)` => `(c₁×c₂)×x`  for `c₁,c₂` constants

These rules are generally safe in the sense of obeying IEEE floating point arithmetic rules. However if the runtime value of `x` happens to be NaN or Inf the **FD** expression `x*0` will identically return 0, because it will have been rewritten to 0 by the simplification rules. The expected IEEE result is NaN.



[^1]: More rules may be added in future versions of FD to improve efficiency.
[^2]: See the [D* ](https://www.microsoft.com/en-us/research/publication/the-d-symbolic-differentiation-algorithm/) paper for an explanation of derivative graph factorization. 
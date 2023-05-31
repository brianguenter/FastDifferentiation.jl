using FastSymbolicDifferentiation

@variables t

A = [t t^2; 3t^2 5]

derivative(A, t)

derivative(A, t, t)


"""
```
julia> A = [t t^2; 3t^2 5]
2×2 Matrix{Node}:
 t              (t ^ 2)
 (3 * (t ^ 2))  5

julia> 

julia> derivative(A, t)
2×2 Matrix{Node}:
 1.0      (2 * t)
 (6 * t)  0.0

julia> 

julia> derivative(A, t, t)
2×2 Matrix{Node{T, 0} where T}:
 0.0  2
 6    0.0

julia> 
```
"""
result() = nothing #this function is here just so you can use tooltips in VSCode to see nicely formatted output  from executing the example code

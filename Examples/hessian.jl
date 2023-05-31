using FastSymbolicDifferentiation

@variables x y z

hessian(x^2 + y^2 + z^2, [x, y, z])

hessian(x * y * z, [x, y, z])

"""
```
julia> x,y,z = Node.((x,y,z))
(x, y, z)

julia> 

julia> hessian(x^2+y^2+z^2,[x,y,z])
3×3 Matrix{Node}:
 2    0.0  0.0
 0.0  2    0.0
 0.0  0.0  2

 julia> hessian(x*y*z,[x,y,z])
 3×3 Matrix{Node}:
  0.0  z    y
  z    0.0  x
  y    x    0.0
```
"""
result() = nothing #this function is here just so you can use tooltips in VSCode to see nicely formatted output  from executing the example code

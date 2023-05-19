using FastSymbolicDifferentiation
using Symbolics

@variables x y z

nx, ny, nz = Node.((x, y, z))

hessian(nx^2 + ny^2 + nz^2, [nx, ny, nz])

hessian(nx * ny * nz, [nx, ny, nz])

"""
```
julia> nx,ny,nz = Node.((x,y,z))
(x, y, z)

julia> 

julia> hessian(nx^2+ny^2+nz^2,[nx,ny,nz])
3×3 Matrix{Node}:
 2    0.0  0.0
 0.0  2    0.0
 0.0  0.0  2

 julia> hessian(nx*ny*nz,[nx,ny,nz])
 3×3 Matrix{Node}:
  0.0  z    y
  z    0.0  x
  y    x    0.0
```
"""
result() = nothing #this function is here just so you can use tooltips in VSCode to see nicely formatted output  from executing the example code

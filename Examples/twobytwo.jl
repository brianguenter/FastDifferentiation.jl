using FastSymbolicDifferentiation
using Symbolics

@variables x y

nx = Node(x)
ny = Node(y)
f1 = cos(nx) * ny
f2 = sin(ny) * nx

gr = DerivativeGraph([f1, f2])
symb = symbolic_jacobian(gr) #non-destructive
func = jacobian_function(gr, [nx, ny])

func(1.0, 2.0)

"""
```
julia> nx = Node(x)
x

julia> ny = Node(y)
y

julia> f1 = cos(nx) * ny
(cos(x) * y)

julia> f2 = sin(ny) * nx
(sin(y) * x)

julia> 

julia> gr = DerivativeGraph([f1, f2])
DerivativeGraph{Int64}(Dict{Node, Int64}((cos(x) * y) => 4, (sin(y) * x) => 6, y => 3, x => 1, sin(y) => 5, cos(x) => 2), Node[x, cos(x), y, (cos(x) * y), sin(y), (sin(y) * x)], Node[(cos(x) * y), (sin(y) * x)], Node[x, y], [4, 6], Dict(4 => 1, 6 => 2), [1, 3], Dict(3 => 2, 1 => 1), Dict{Int64, FastSymbolicDifferentiation.EdgeRelations{Int64}}(5 => FastSymbolicDifferentiation.EdgeRelations{Int64}(FastSymbolicDifferentiation.PathEdge{Int64}[(6 5  1 x Bool[0, 1] Bool[0, 1])], FastSymbolicDifferentiation.PathEdge{Int64}[(5 3  1 cos(y) Bool[0, 1] Bool[0, 1])]), 4 => FastSymbolicDifferentiation.EdgeRelations{Int64}(FastSymbolicDifferentiation.PathEdge{Int64}[], FastSymbolicDifferentiation.PathEdge{Int64}[(4 2  1 y Bool[1, 0] Bool[1, 0]), (4 3  1 cos(x) Bool[1, 0] Bool[0, 1])]), 6 => FastSymbolicDifferentiation.EdgeRelations{Int64}(FastSymbolicDifferentiation.PathEdge{Int64}[], FastSymbolicDifferentiation.PathEdge{Int64}[(6 5  1 x Bool[0, 1] Bool[0, 1]), (6 1  1 sin(y) Bool[0, 1] Bool[1, 0])]), 2 => FastSymbolicDifferentiation.EdgeRelations{Int64}(FastSymbolicDifferentiation.PathEdge{Int64}[(4 2  1 y Bool[1, 0] Bool[1, 0])], FastSymbolicDifferentiation.PathEdge{Int64}[(2 1  1 -(sin(x)) Bool[1, 0] Bool[1, 0])]), 3 => FastSymbolicDifferentiation.EdgeRelations{Int64}(FastSymbolicDifferentiation.PathEdge{Int64}[(4 3  1 cos(x) Bool[1, 0] Bool[0, 1]), (5 3  1 cos(y) Bool[0, 1] Bool[0, 1])], FastSymbolicDifferentiation.PathEdge{Int64}[]), 1 => FastSymbolicDifferentiation.EdgeRelations{Int64}(FastSymbolicDifferentiation.PathEdge{Int64}[(2 1  1 -(sin(x)) Bool[1, 0] Bool[1, 0]), (6 1  1 sin(y) Bool[0, 1] Bool[1, 0])], FastSymbolicDifferentiation.PathEdge{Int64}[])), IdDict{Any, Any}())

julia> symb = symbolic_jacobian(gr) #non-destructive
2×2 Matrix{Node}:
 (y * -(sin(x)))  cos(x)
 sin(y)           (x * cos(y))

julia> func = jacobian_function(gr, [nx, ny])
RuntimeGeneratedFunction(#=in FastSymbolicDifferentiation=#, #=using FastSymbolicDifferentiation=#, :((x, y)->begin
          result = fill(0.0, (2, 2))
          begin
              var"##295" = sin(x)
              var"##294" = -var"##295"
              var"##293" = y * var"##294"
              result[CartesianIndex(1, 1)] = var"##293"
          end
          begin
              var"##296" = sin(y)
              result[CartesianIndex(2, 1)] = var"##296"
          end
          begin
              var"##297" = cos(x)
              result[CartesianIndex(1, 2)] = var"##297"
          end
          begin
              var"##299" = cos(y)
              var"##298" = x * var"##299"
              result[CartesianIndex(2, 2)] = var"##298"
          end
          return result
      end))

julia> 

julia> func(1.0, 2.0)
2×2 Matrix{Float64}:
 -1.68294    0.540302
  0.909297  -0.416147
```
"""
result() = nothing #this function is here just so you can use tooltips in VSCode to see nicely formatted output  from executing the example code

using FastSymbolicDifferentiation

@variables x y

f1 = cos(x) * y
f2 = sin(y) * x

symb = jacobian([f1, f2], [x, y]) #non-destructive
func = jacobian_function([f1, f2], [x, y])

func(1.0, 2.0)

"""
```
julia> @variables x y

julia> f1 = cos(x) * y
(cos(x) * y)

julia> f2 = sin(y) * x
(sin(y) * x)

julia> 

julia> symb = jacobian([f1, f2], [x, y]) #non-destructive
2×2 Matrix{Node}:
 (y * -(sin(x)))  cos(x)
 sin(y)           (x * cos(y))

julia> func = jacobian_function([f1, f2], [x, y])
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
result() = nothing #this function is here just so you can use tooltips in VSCode to see nicely formatted output from executing the example code

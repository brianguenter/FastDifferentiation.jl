using FastSymbolicDifferentiation

@variables x y

f1 = cos(x) * y
f2 = sin(y) * x

symb = jacobian([f1, f2], [x, y]) #non-destructive
func = make_function(symb, [x, y])

func([1.0, 2.0])

"""
```
julia> @variables x y
y

julia> f1 = cos(x) * y
(cos(x) * y)

julia> f2 = sin(y) * x
(sin(y) * x)

julia> symb = jacobian([f1, f2], [x, y]) #non-destructive
2×2 Matrix{Node}:
 (y * -(sin(x)))  cos(x)
 sin(y)           (x * cos(y))

julia> func = make_function(symb, [x, y])
RuntimeGeneratedFunction(#=in FastSymbolicDifferentiation=#, #=using FastSymbolicDifferentiation=#, :((input_variables,)->begin
          #= /home/brian/source/FastSymbolicDifferentiation.jl/src/CodeGeneration.jl:25 =#
          begin
              result = Array{promote_type(Float64, eltype(input_variables)), 2}(undef, 2, 2)
              begin
                  var"##295" = sin(input_variables[1])
                  var"##294" = -var"##295"
                  var"##293" = input_variables[2] * var"##294"
                  result[CartesianIndex(1, 1)] = var"##293"
              end
              begin
                  var"##296" = sin(input_variables[2])
                  result[CartesianIndex(2, 1)] = var"##296"
              end
              begin
                  var"##297" = cos(input_variables[1])
                  result[CartesianIndex(1, 2)] = var"##297"
              end
              begin
                  var"##299" = cos(input_variables[2])
                  var"##298" = input_variables[1] * var"##299"
                  result[CartesianIndex(2, 2)] = var"##298"
              end
              return result
          end
      end))

julia> func([1.0, 2.0])
2×2 Matrix{Float64}:
 -1.68294    0.540302
  0.909297  -0.416147
```
"""
result() = nothing #this function is here just so you can use tooltips in VSCode to see nicely formatted output from executing the example code

# How to use `make_function``
## `make_function` options

The `make_function` procedure creates a runtime generated function from **FD** expressions. The signature for `make_function` is

`make_function(func_array::AbstractArray{T}, input_variables::AbstractVector{<:Node}...; in_place::Bool=false, init_with_zeros::Bool=true) where {T<:Node}`

The type of the first argument, `func_array`, determines the type of the array that will hold the result of evaluating the **FD** expressions. These are rough rules of thumb for parameter settings for maximum performance:
* If `func_array` is dense and has more than 100 elements then `typeof(func_array)` should be a subtype of `Array` with `in_place=true`, `init_with_zeros=false` (but see caveat below)
* If `func_array` is sparse (use the [`sparsity`](@ref) function to measure this) then `typeof(func_array)` should be `SparseMatrixCSC` (automatically created by [`sparse_jacobian`](@ref) and [`sparse_hessian`](@ref) functions) and `in_place=true`. The `init_with_zeros` argument has no effect for sparse arrays.
* if `func_array` is small, less than 100 elements, then `typeof(func_array)` should be a subtype of StaticArray, `in_place=false`, `init_with_zeros=true`.

The type of the second argument `input_variables` has no effect on performance of the runtime generated function.

If the third argument, `in_place`, = false then the generated code will create and return a zero initialized array at every call. If `in_place = true` the generated code has two arguments: an array to accept the result of the evaluation, and a vector of inputs. There will only be a performance advantage by setting `in_place = true` if the type of `func_array` is `Array` or `SparseMatrixCSC`.

The fourth argument, `init_with_zeros`, only has an effect when `in_place=true` *and* the type of `func_array` is `Array` or `SArray`. If `init_with_zeros = true` then the in place array argument will be initialized with zeros at each call. If it is false it will not be initialized with zeros. 

This is useful when you are evaluating a function many times but not modifying the result array. You can initialize the in place array with zeros once before passing it to the runtime generated function. The runtime generated function will not alter identically zero elements so they will stay zero. Be careful though that some other function doesn't modify the in place array between calls to the runtime generated function. This could corrupt the identically zero elements.

Example: generated function returns static array
```julia

@variables x y

julia> func = SA[x^2+y,y*x,y^2]
3-element SVector{3, FastDifferentiation.Node} with indices SOneTo(3):
 ((x ^ 2) + y)
       (y * x)
       (y ^ 2)

julia> f_exe! = make_function(func,[x,y])
...

julia> f_exe!([2.0,3.0])
3-element SVector{3, Float64} with indices SOneTo(3):
 7.0
 6.0
 9.0
```

Example: assume your result is a large dense array (> 100 elements) and that you are using an in place array with no initialization. For dense arrays this should generate the fastest code:
```
julia> @variables x y z
z

julia> p = [x^2 y^2; z^2 0]
2×2 Matrix{FastDifferentiation.Node}:
 (x ^ 2)  (y ^ 2)
 (z ^ 2)        0

julia> floatp = similar(p,Float64)
2×2 Matrix{Float64}:
 5.0e-324  1.08279e-311
 1.5e-323  1.08279e-311

julia> fexe! = make_function(p,[x,y,z],in_place=true,init_with_zeros=false)
...


julia> fexe!(floatp,[1.0,2.0,3.0])
4.0

julia> floatp
2×2 Matrix{Float64}:
 1.0  4.0
 9.0  1.08279e-311

```

Notice that `floatp[2,2]` is not 0. This is because `init_with_zeros=false`. If you need this array element to be zero then initialize it before making the call to `fexe!`.


You can use the [`sparsity`](@ref) function to measure sparsity and determine whether you should use the dense or sparse derivative functions.

## Evaluate a function and derivatives
Sometimes you want to evaluate a function and one or more derivative orders. If you pack all the terms you want to evaluate into the argument to `make_function` then common terms will be detected and only computed once. This will be generally be more efficient than evaluating the function and derivatives separately:

```julia
julia> f = [x^2*y^2,sqrt(x*y)]
2-element Vector{FastDifferentiation.Node}:
 ((x ^ 2) * (y ^ 2))
       sqrt((x * y))

julia> jac = jacobian(f,[x,y])
2×2 Matrix{FastDifferentiation.Node}:
             ((y ^ 2) * (2 * x))              ((x ^ 2) * (2 * y))
 ((1 / (2 * sqrt((x * y)))) * y)  ((1 / (2 * sqrt((x * y)))) * x)


julia> f_and_jac = make_function([vec(jac);f],[x,y])
...

julia> tmp = f_and_jac([1.1,2.1])
6-element Vector{Float64}:
 9.702000000000002
 0.6908492797077573
 5.082000000000001
 0.36187343222787294
 5.336100000000001
 1.5198684153570665

julia> jac_eval = reshape(view(tmp,1:4),2,2)
2×2 reshape(view(::Vector{Float64}, 1:4), 2, 2) with eltype Float64:
 9.702     5.082
 0.690849  0.361873

julia> f_eval = view(tmp,5:6)
2-element view(::Vector{Float64}, 5:6) with eltype Float64:
 5.336100000000001
 1.5198684153570665
```

There are several options for `make_function`. If `in_place==false`, the default, then it will create and return a new matrix at each function call. If `in_place==true` it will make a function that expects two arguments, a matrix to hold the result and a vector of input variable values. The `in_place` option is available on all executables including Jᵀv,Jv,Hv.

```julia
julia> jac_exe! = make_function(symb,[x,y], in_place=true)
...
julia> a = similar(symb,Float64)
2×2 Matrix{Float64}:
 0.0  0.0
 0.0  6.93532e-310

julia> jac_exe!(a,[1.0,2.0])
2×2 Matrix{Float64}:
 -1.68294    0.540302
  0.909297  -0.416147

julia> a
2×2 Matrix{Float64}:
 -1.68294    0.540302
  0.909297  -0.416147
```

## Generating code
If you need to generate code that can be cut and pasted into another application you can use `make_Expr` instead of `make_function`. It has the same arguments except the `in_place` and `init_with_zeros` boolean arguments are not optional.
```julia
julia> @variables x y
y

julia> j = jacobian([x^2*y^2,cos(x+y),log(x/y)],[x,y])
3×2 Matrix{FastDifferentiation.Node}:
       ((y ^ 2) * (2 * x))                 ((x ^ 2) * (2 * y))
           -(sin((x + y)))                     -(sin((x + y)))
 ((1 / (x / y)) * (1 / y))  ((1 / (x / y)) * -(((x / y) / y)))

julia> in_place_no_init = make_Expr(j,[x,y],true,false)
:((result, input_variables)->begin
          #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:127 =#
          #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:127 =# @inbounds begin
                  #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:128 =#
                  begin
                      var"##431" = input_variables[2] ^ 2
                      var"##432" = 2 * input_variables[1]
                      var"##430" = var"##431" * var"##432"
                      result[CartesianIndex(1, 1)] = var"##430"
                      var"##435" = input_variables[1] + input_variables[2]
                      var"##434" = sin(var"##435")
                      var"##433" = -var"##434"
                      result[CartesianIndex(2, 1)] = var"##433"
                      var"##438" = input_variables[1] / input_variables[2]
                      var"##437" = 1 / var"##438"
                      var"##439" = 1 / input_variables[2]
                      var"##436" = var"##437" * var"##439"
                      result[CartesianIndex(3, 1)] = var"##436"
                      var"##441" = input_variables[1] ^ 2
                      var"##442" = 2 * input_variables[2]
                      var"##440" = var"##441" * var"##442"
                      result[CartesianIndex(1, 2)] = var"##440"
                      result[CartesianIndex(2, 2)] = var"##433"
                      var"##445" = var"##438" / input_variables[2]
                      var"##444" = -var"##445"
                      var"##443" = var"##437" * var"##444"
                      result[CartesianIndex(3, 2)] = var"##443"
                  end
              end
      end)

julia> in_place_zero_init = make_Expr(j,[x,y],true,true)
:((result, input_variables)->begin
          #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:127 =#
          #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:127 =# @inbounds begin
                  #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:128 =#
                  begin
                      var"##447" = input_variables[2] ^ 2
                      var"##448" = 2 * input_variables[1]
                      var"##446" = var"##447" * var"##448"
                      result[CartesianIndex(1, 1)] = var"##446"
                      var"##451" = input_variables[1] + input_variables[2]
                      var"##450" = sin(var"##451")
                      var"##449" = -var"##450"
                      result[CartesianIndex(2, 1)] = var"##449"
                      var"##454" = input_variables[1] / input_variables[2]
                      var"##453" = 1 / var"##454"
                      var"##455" = 1 / input_variables[2]
                      var"##452" = var"##453" * var"##455"
                      result[CartesianIndex(3, 1)] = var"##452"
                      var"##457" = input_variables[1] ^ 2
                      var"##458" = 2 * input_variables[2]
                      var"##456" = var"##457" * var"##458"
                      result[CartesianIndex(1, 2)] = var"##456"
                      result[CartesianIndex(2, 2)] = var"##449"
                      var"##461" = var"##454" / input_variables[2]
                      var"##460" = -var"##461"
                      var"##459" = var"##453" * var"##460"
                      result[CartesianIndex(3, 2)] = var"##459"
                  end
              end
      end)
```
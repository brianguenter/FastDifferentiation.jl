# Examples


**FD** uses a global cache for common subexpression elimination so the **FD** symbolics preprocessing step is not thread safe. 

Under ordinary conditions the memory used by the cache won't be an issue. But, if you have a long session where you are creating many complex functions it is possible the cache will use too much memory. If this happens call the function `clear_cache` after you have completely processed your expression.

The most common way to use **FD** is this:
* create variables
* do operations on those variables to create the function you want to differentiate
* compute a symbolic derivative of the function
* pass the symbolic derivative to `make_function` to generate a function to efficiently evaluate the derivative


##### Creating Variables

Scalar variables
```julia
using FastDifferentiation

@variables x y z

```
Arrays of variables of arbitrary dimension
```julia
julia> make_variables(:x,3)
3-element Vector{FastDifferentiation.Node}:
 x1
 x2
 x3

julia> make_variables(:x,2,3)
2×3 Matrix{FastDifferentiation.Node}:
 x1_1  x1_2  x1_3
 x2_1  x2_2  x2_3

julia> make_variables(:x,2,3,2)
2×3×2 Array{FastDifferentiation.Node, 3}:
[:, :, 1] =
 x1_1_1  x1_2_1  x1_3_1
 x2_1_1  x2_2_1  x2_3_1

[:, :, 2] =
 x1_1_2  x1_2_2  x1_3_2
 x2_1_2  x2_2_2  x2_3_2
```

##### Compute derivatives

Compute higher order derivatives
```julia
julia> @variables x y
y

julia> f = x^3*y^3
((x ^ 3) * (y ^ 3))

julia> derivative([f],x,y,x) #take derivative wrt x, then y, then x
1-element Vector{FastDifferentiation.Node{typeof(*), 2}}:
 (18 * (x * (y ^ 2)))

 julia> derivative([cos(x*y);;;exp(x*y)],x,y,x) #derivative accepts input arrays of any dimension
1×1×2 Array{FastDifferentiation.Node{typeof(+), 2}, 3}:
[:, :, 1] =
 ((-(y) * cos((x * y))) + ((((x * -(y)) * -(sin((x * y)))) + -(cos((x * y)))) * y))

[:, :, 2] =
 (((((x * y) + 1) * exp((x * y))) * y) + (y * exp((x * y))))
```

Compute derivative of a function and make executable

```julia
# compute Jacobian and generate function to evaluate it
julia> f1 = cos(x) * y
(cos(x) * y)

julia> f2 = sin(y) * x
(sin(y) * x)

julia> symb = jacobian([f1, f2], [x, y]) #the vector [x,y] tells make_function 
# how to order the arguments to the generated function
2×2 Matrix{Node}:
 (y * -(sin(x)))  cos(x)
 sin(y)           (x * cos(y))

julia> jac_exe = make_function(symb,[x,y]) 
...
julia> jac_exe([1.0,2.0]) #jac_exe was created with variable ordering [x,y] 
# so x will get the value 1.0, y 2.0
2×2 Matrix{Float64}:
 -1.68294    0.540302
  0.909297  -0.416147
```

```julia
# compute Hessian and generate function to evaluate it
@variables x y z

julia> h_symb1 = hessian(x^2*y^2*z^2,[x,y,z])
3×3 Matrix{FastDifferentiation.Node}:
 (2 * ((z ^ 2) * (y ^ 2)))        (((2 * x) * (2 * y)) * (z ^ 2))  (((2 * x) * (2 * z)) * (y ^ 2))
 (((2 * y) * (2 * x)) * (z ^ 2))  (2 * ((z ^ 2) * (x ^ 2)))        (((2 * y) * (2 * z)) * (x ^ 2))
 (((2 * z) * (2 * x)) * (y ^ 2))  (((2 * z) * (2 * y)) * (x ^ 2))  (2 * ((x ^ 2) * (y ^ 2)))

julia> hexe_1 = make_function(h_symb1,[x,y,z]) #the vector [x,y,z] tells make_function 
# how to order the arguments to the generated function
...
julia> hexe_1([1.0,2.0,3.0]) #hexe_1 was created with variable ordering [x,y,z] 
# so x will get the value 1.0, y 2.0, and z 3.0
3×3 Matrix{Float64}:
 72.0  72.0  48.0
 72.0  18.0  24.0
 48.0  24.0   8.0
```


Compute any subset of the columns of the Jacobian:
```julia
julia> symb = jacobian([x*y,y*z,x*z],[x,y,z]) #all columns
3×3 Matrix{Node}:
 y    x    0.0
 0.0  z    y
 z    0.0  x

julia> symb = jacobian([x*y,y*z,x*z],[x,y]) #first two columns
3×2 Matrix{Node}:
 y    x
 0.0  z
 z    0.0

julia> symb = jacobian([x*y,y*z,x*z],[z,y]) #second and third columns, ordered so ∂f/∂z is 1st column of the output, ∂f/∂y the 2nd
3×2 Matrix{Node}:
 0.0  x
 y    z
 x    0.0
```
##### More on make_function
There are many options to `make_function`. This brief summary will help you figure out how to set them to achieve maximum performance.

There are two categories of problem for which different settings are likely to yield best performance: small and large. The first category has a return matrix that is small with less than 100 elements, and an input vector that in this size range. In this case you will probably get the fastest results by calling `make_function` with a static array as the first argument and then calling the runtime generated function with an SVector or MVector for the input variable arguments.

If your problem is in the large category then whenever possible set `in_place=true` and pass a matrix into the runtime generated function to hold the return values. You can improve performance further by combining `in_place=true` with `init_with_zeros=false`. The runtime generated function will not initialize your in place array with zeros each time it is called. 

You can use the [`sparsity`](@ref) function to determine if any entries in your function are identically zero. If none are then `init_with_zeros=false` is safe even if you don't pre-initialize your array with zeros because every array entry will be set with a valid value.

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

For out of place evaluation (matrix created and returned by the executable function) the input vector and return matrix of the executable can be any mix of StaticArray and Vector. If the first argument to `make_function` is a subtype of StaticArray then the compiled executable will return a StaticArray value. The compiled executable can be called with either an `SVector` or `Vector` argument. For small input sizes the `SVector` should be faster, essentially the same as passing the input as scalar values.

For functions with low input and output dimensions the fastest executable will be generated by calling `make_function` with first argument a subtype of StaticArray and calling the executable with an SVector argument. The usual cautions of StaticArrays apply, that total length of the return value < 100 or so and total length of the input < 100 or so.

```julia
julia> @variables x y
y

julia> j = jacobian([x^2 * y^2, cos(x + y), log(x / y)], [x, y])

julia> j_exe = make_function(j, [x, y])

julia> @assert typeof(j_exe([1.0, 2.0])) <: Array #return type is Array and input type is Vector

julia> j_exe2 = make_function(SArray{Tuple{3,2}}(j), [x, y])

julia> @assert typeof(j_exe2(SVector{2}([1.0, 2.0]))) <: StaticArray #return type is StaticArray and input type is SVector. This should be the fastest.
```
If you need to generate code that can be cut and pasted into another application then you can use `make_Expr` instead of `make_function`. It has the same arguments except the `in_place` and `init_with_zeros` boolean arguments are not optional.
```julia
julia> @variables x y
y

julia> j = jacobian([x^2 * y^2, cos(x + y), log(x / y)], [x, y])

3×2 Matrix{FastDifferentiation.Node}:
 ((y ^ 2) * (2 * x))           ((x ^ 2) * (2 * y))
     -(sin((x + y)))               -(sin((x + y)))
 ((y / x) * (1 / y))  ((y / x) * -(((x / y) / y)))

#create code to evaluate the jacobian with in place matrix
julia> make_Expr(j,[x,y],true,false) 
#Notice that the generated code does NOT initialize the input to 0
#
:((result, input_variables)->begin
          #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:94 =#
          #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:94 =# @inbounds begin
                  #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:95 =#
                  begin
                      begin
                          var"##293" = input_variables[2] ^ 2
                          var"##294" = 2 * input_variables[1]
                          var"##292" = var"##293" * var"##294"
                          result[CartesianIndex(1, 1)] = var"##292"
                      end
                      begin
                          var"##297" = input_variables[1] + input_variables[2]
                          var"##296" = sin(var"##297")
                          var"##295" = -var"##296"
                          result[CartesianIndex(2, 1)] = var"##295"       
                      end
                      begin
                          var"##300" = input_variables[1] / input_variables[2]
                          var"##299" = 1 / var"##300"
                          var"##301" = 1 / input_variables[2]
                          var"##298" = var"##299" * var"##301"
                          result[CartesianIndex(3, 1)] = var"##298"       
                      end
                      begin
                          var"##303" = input_variables[1] ^ 2
                          var"##304" = 2 * input_variables[2]
                          var"##302" = var"##303" * var"##304"
                          result[CartesianIndex(1, 2)] = var"##302"       
                      end
                      begin
                          result[CartesianIndex(2, 2)] = var"##295"       
                      end
                      begin
                          var"##307" = var"##300" / input_variables[2]    
                          var"##306" = -var"##307"
                          var"##305" = var"##299" * var"##306"
                          result[CartesianIndex(3, 2)] = var"##305"       
                      end
                  end
              end
      end)

julia> make_Expr(j, [x, y], true, true) 
#Notice that the generated code initializes the input array to zero
:((result, input_variables)->begin
          #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:94 =#
          #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:94 =# @inbounds begin
                  #= c:\Users\seatt\source\FastDifferentiation.jl\src\CodeGeneration.jl:95 =#
                  begin
                      result .= zero(eltype(input_variables)) #INITIALIZE to 0
                      begin
                          var"##309" = input_variables[2] ^ 2
                          var"##310" = 2 * input_variables[1]
                          var"##308" = var"##309" * var"##310"
                          result[CartesianIndex(1, 1)] = var"##308"       
                      end
                      begin
                          var"##313" = input_variables[1] + input_variables[2]
                          var"##312" = sin(var"##313")
                          var"##311" = -var"##312"
                          result[CartesianIndex(2, 1)] = var"##311"       
                      end
                      begin
                          var"##316" = input_variables[1] / input_variables[2]
                          var"##315" = 1 / var"##316"
                          var"##317" = 1 / input_variables[2]
                          var"##314" = var"##315" * var"##317"
                          result[CartesianIndex(3, 1)] = var"##314"       
                      end
                      begin
                          var"##319" = input_variables[1] ^ 2
                          var"##320" = 2 * input_variables[2]
                          var"##318" = var"##319" * var"##320"
                          result[CartesianIndex(1, 2)] = var"##318"       
                      end
                      begin
                          result[CartesianIndex(2, 2)] = var"##311"       
                      end
                      begin
                          var"##323" = var"##316" / input_variables[2]    
                          var"##322" = -var"##323"
                          var"##321" = var"##315" * var"##322"
                          result[CartesianIndex(3, 2)] = var"##321"       
                      end
                  end
              end
      end)
```
##### Sparse Jacobians and Hessians
The functions `sparse_jacobian, sparse_hessian` compute sparse symbolic derivatives. When you pass a sparse symbolic function matrix to `make_function` it will generate an executable which expects an in place sparse matrix to hold the result. For functions with sparse Jacobians or Hessians this can be orders of magnitude faster than using a dense in place matrix.

```julia
julia> hess = sparse_hessian(x^3 + y^3 + z^3, [x,y,z])        
3×3 SparseArrays.SparseMatrixCSC{FastDifferentiation.Node, Int64} with 3 stored entries:
 (6 * x)        ⋅        ⋅
       ⋅  (6 * y)        ⋅
       ⋅        ⋅  (6 * z)


julia> res = similar(hess,Float64) #make sparse matrix with proper sparsity to pass to the generated function
3×3 SparseArrays.SparseMatrixCSC{Float64, Int64} with 3 stored entries:
 0.0   ⋅    ⋅ 
  ⋅   0.0   ⋅
  ⋅    ⋅   0.0

julia> sp_f = make_function(hess,[x,y,z])
...

julia> sp_f([1.0,2.0,3.0],res)
3×3 SparseArrays.SparseMatrixCSC{Float64, Int64} with 3 stored entries:
 6.0    ⋅     ⋅
  ⋅   12.0    ⋅
  ⋅     ⋅   18.0
```
##### Less commonly used functions
Compute `Hv` without forming the full Hessian matrix. This is useful if the Hessian is very large
```julia
julia> @variables x y
y

julia> f = x^2 * y^2
((x ^ 2) * (y ^ 2))

julia> hv_fast, v_vec2 = hessian_times_v(f, [x, y])
...

julia> hv_fast_exe = make_function(hv_fast, [[x, y]; v_vec2]) #need v_vec2 because hv_fast is a function of x,y,v1,v2 and have to specify the order of all inputs to the executable
...
julia> hv_fast_exe([1.0,2.0,3.0,4.0]) #first two vector elements are x,y last two are v1,v2
2-element Vector{Float64}:
 56.0
 32.0
```

Symbolic and executable Jᵀv and Jv (see this [paper](https://arxiv.org/abs/1812.01892) for applications of this operation).
```julia
julia> (f1,f2) = cos(x)*y,sin(y)*x
((cos(x) * y), (sin(y) * x))

julia> jv,vvec = jacobian_times_v([f1,f2],[x,y])
...

julia> jv_exe = make_function(jv,[[x,y];vvec])
...

julia> jv_exe([1.0,2.0,3.0,4.0]) #first 2 arguments are x,y values and last two are v vector values

2×1 Matrix{Float64}:
 -2.8876166853748195
  1.0633049342884753

julia> jTv,rvec = jacobian_transpose_v([f1,f2],[x,y])
...

julia> jtv_exe = make_function(jTv,[[x,y];rvec])
...
julia> jtv_exe([1.0,2.0,3.0,4.0])
2-element Vector{Float64}:
 -1.4116362015446517
 -0.04368042858415033
```
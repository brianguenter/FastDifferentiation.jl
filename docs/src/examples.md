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
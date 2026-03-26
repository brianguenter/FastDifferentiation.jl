# API

## Lazy Differentiation

Create unevaluated derivative nodes with [`lazy_differential`](@ref), then resolve them with [`expand_derivatives`](@ref) before passing to solvers:

```julia
@variables x y
f = x^2 * y

Dt = lazy_differential(x)
lazy_expr = Dt(f)                        # unevaluated ∂f/∂x
expanded  = expand_derivatives(lazy_expr) # computes 2x*y

jac = jacobian([expanded], [x, y])
```

Nested and multi-variable lazy derivatives are supported:

```julia
dx = lazy_differential(x)
dy = lazy_differential(y)
nested = dy(dx(f))                        # ∂²f/∂y∂x (unevaluated)
expanded = expand_derivatives(nested)     # computes 2x
```

```@docs
lazy_differential
expand_derivatives
```

## All Exported Functions

```@autodocs
Modules = [FastDifferentiation]
Private = false
Order = [:macro, :function]
```
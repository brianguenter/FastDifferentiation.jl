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

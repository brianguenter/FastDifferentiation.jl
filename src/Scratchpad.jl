#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables x y

    uf = UnspecifiedFunction(:q, SVector(Node(x), Node(y)), SVector{0,Node{SymbolicUtils.BasicSymbolic{Real},0}}())
    ufn = Node(uf, uf.variables...)
    deriv = node_value(derivative(derivative(ufn, Val{1}()), Val{2}()))

    @assert deriv.variables == SVector(Node(x), Node(y))

    fn = DerivativeGraph(x * ufn)
    jac = symbolic_jacobian!(fn)
    println(jac[1, 1])
    @assert jac[1, 1] == ufn + x * derivative(ufn, Val{1}())
    @assert jac[1, 2] == x * derivative(ufn, Val{2}())

end
export test

#changed
#change
export test

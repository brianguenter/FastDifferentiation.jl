#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables q1 q2
    nq1 = Node(q1)
    nq2 = Node(q2)

    A = [
        cos(nq1) -cos(nq1)
        sin(nq1) sin(nq1)
    ]

    DA = Node[
        -sin(nq1) sin(nq1)
        cos(nq1) cos(nq1)
    ]

    @assert isapprox(zeros(2, 2), node_value.(derivative(A, nq2))) #taking derivative wrt variable not present in the graph returns all zero matrix
    println(derivative(A, nq1))
    println(DA)
    @assert DA == derivative(A, nq1)




end
export test

#changed
#change
export test

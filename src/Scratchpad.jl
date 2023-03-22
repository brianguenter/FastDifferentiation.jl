#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables q1, q2
    nq1 = Node(q1)
    nq2 = Node(q2)

    A = [
        cos(nq1) -cos(nq1)
        sin(nq1) sin(nq1)
    ]

    DerivativeGraph([-cos(nq1)])
    # derivative(A, q1, q1)
    derivative(A, nq1)
end
export test

#changed
#change
export test

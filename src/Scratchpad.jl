#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables x, y
    uf = UnspecifiedFunction(:q, SVector(Node(x), Node(y)), SVector(Node(x)))



end
export test

#changed
#change
export test

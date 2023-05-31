#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences



function test()
    @variables x y


    a = x * y
    @assert derivative(a, Val(1)) == y
    @assert derivative(a, Val(2)) == x
    @assert derivative(x) == Node(1)
    @assert derivative(Node(1)) == Node(0)

end




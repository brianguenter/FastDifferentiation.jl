#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    @variables x
    nx = Node(x)
    tmp0 = make_function([nx * nx], [nx])
    # dfsimp(x) = tmp0([x])[1]
end




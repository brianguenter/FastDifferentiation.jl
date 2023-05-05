#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    _, graph, _, _ = simple_dominator_graph()
    factor!(graph)
end
export test



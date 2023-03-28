#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests

function test()
    @variables x y



    q = function_of(:q, x, y)
    f = x * q + y * q
    graph = DerivativeGraph([f])
    # subs, subs_dict = compute_factorable_subgraphs(graph)
    # Vis.draw(graph)
    # readline()
    # for sub in subs
    #     println("factored subgraph $sub")
    #     factor_one_subgraph!(graph, sub, subs, subs_dict)
    #     Vis.draw(graph)
    #     readline()
    # end
    symbolic_jacobian!(graph)
end
export test

#changed
#change
export test

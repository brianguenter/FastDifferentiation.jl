#this file is for temporary testing code since it is so hard to debug tests using the VSCode test system. 

using Symbolics
using StaticArrays
using FiniteDifferences
using .FSDTests


function test()
    @variables x, y, z

    fsd_graph = spherical_harmonics(FastSymbolic(), 10, x, y, z)
    sprse = sparse_symbolic_jacobian!(fsd_graph, variables(fsd_graph))
    fsd_graph = spherical_harmonics(FastSymbolic(), 5, x, y, z) #because global cache has not been reset the sparse and dense graphs should have identical elements.
    dense = symbolic_jacobian!(fsd_graph, variables(fsd_graph))

    # for index in CartesianIndices(sprse)
    #     @test sprse[index] == dense[index]
    # end

    for index in CartesianIndices(dense)
        if sprse[index] != dense[index] #empty elements in sprse get value Node{Int64,0} wherease zero elements in dense get value Node{Float64,0}. These are not == so need special case.
            @assert value(sprse[index]) == value(dense[index])
        else
            @assert sprse[index] == dense[index]
        end
    end
end
export test



module NonUnique

using JLD2: JLD2
using FastDifferentiation
import FastDifferentiation as FD
# using ElectronDisplay


function bug()
    erroneous_inputs = JLD2.load_object("fd_error_inputs.jld2")

    (; f_node, z_node) = erroneous_inputs

    # root_node = 81
    gr = FD.DerivativeGraph(f_node)
    # edges = FD.child_edges(gr, 40)
    # nd = FD.node(gr, root_node)
    # subgr = FD.DerivativeGraph([nd])




    FD.factor!(gr)
    edg = FD.edges(gr, 87, 40)
    display(edg)
    # FD.factor!(subgr)
    FD.write_dot("full.svg", gr, value_labels=true, reachability_labels=false)

    # FD.write_dot("$root_node-bug.svg", subgr, value_labels=false)
    # # this hits an assertion
    # # ERROR: AssertionError: Should only be one path from root 2 to variable 6. Instead have 2 children from node 1474 on the path
    # tmp = FD.DerivativeGraph(f_node)
    # FD.factor!(tmp)
    # FD.write_dot("f_node-bug.svg", tmp, value_labels=false)
    # FD.sparse_jacobian(f_node, z_node)

    return nothing
end
export bug

import FiniteDifferences
using FastDifferentiation.FDInternals
using FastDifferentiation.FDTests
function test()


    FD_graph = spherical_harmonics(6)
    factor!(FD_graph)
    FD.write_dot("sph.svg", FD_graph, value_labels=true, reachability_labels=false)
    # FD.jacobian(roots(FD_graph), variables(FD_graph))
    return nothing
end
export test


end # module NonUnique
export NonUnique


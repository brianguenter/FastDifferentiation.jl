function verify_invariants(v1::Integer, v2::Integer)
    @assert v1 !== v2 "vertex cannot have an edge to itself."
    @assert v1 > 0 && v2 > 0 #nodes must have positive numbers
    #make sure first vertex always has higher number than second.
    if v1 < v2
        tmp = v1
        v1 = v2
        v2 = tmp
    end
    return v1, v2
end

""""
    PathEdge{T<:Integer,N,M}

Represents a partial derivative edge in a derivative graph. Edge keeps track of which roots and variables are reachable from this edge.

# Fields

- `top_vertex::T`
- `bott_vertex::T`
- `edge_value::Node`
- `reachable_variables::BitVector`
- `reachable_roots::BitVector`
- `unique_id::UUID`
"""
struct PathEdge{T<:Integer,N,M}
    top_vertex::T
    bott_vertex::T
    edge_value::Node
    reachable_variables::BitVector
    reachable_roots::BitVector
    unique_id::UUID #This reliably distiguishes between two otherwise identical edges. Otherwise the only difference is that the BitVector arrays of two PathEdges use different physical arrays even if the values are identical. Relying on this to make Set operations work correctly is too subtle and fragile.

    PathEdge(v1::T, v2::T, expression::S, domain_dim::T, codomain_dim::T) where {T<:Integer,S<:Real} = PathEdge(v1, v2, Node(expression), domain_dim, codomain_dim)

    function PathEdge(v1::T, v2::T, expression::Node, domain_dim::T, codomain_dim::T) where {T<:Integer}
        v2, v2 = verify_invariants(v1, v2)

        #Edges need to have unique ID. Otherwise it is not possible to have two vertices connected by more than one edge.
        return new{T,domain_dim,codomain_dim}(v1, v2, expression, BitVector(undef, domain_dim), BitVector(undef, codomain_dim), uuid1())
    end

    PathEdge(v1::T, v2::T, expression::S, reachable_variables::BitVector, reachable_roots::BitVector) where {T,S<:Real} = PathEdge(v1, v2, Node(expression), reachable_variables, reachable_roots, uuid1()) #if expression is a number rather than a Node wrap it in a Node.

    function PathEdge(v1::T, v2::T, expression::Node, reachable_variables::BitVector, reachable_roots::BitVector) where {T}
        N = length(reachable_variables)
        M = length(reachable_roots)
        v1, v2 = verify_invariants(v1, v2)

        return new{T,N,M}(v1, v2, expression, reachable_variables, reachable_roots, uuid1())
    end
end



value(e::PathEdge) = e.edge_value

top_vertex(e::PathEdge) = e.top_vertex

bott_vertex(e::PathEdge) = e.bott_vertex

vertices(e::PathEdge) = (top_vertex(e), bott_vertex(e))


times_used(a::PathEdge) = sum(reachable_roots(a)) * sum(reachable_variables(a))

"""
    reachable_roots(e::PathEdge)

This returns the internal data structure so **DO NOT MODIFY** the returned value unless you know what you are doing! Modification will corrupt the derivative graph. If you need to modify the result then copy first."""
reachable_roots(e::PathEdge) = e.reachable_roots

is_root_reachable(e::PathEdge, root_index::Integer) = reachable_roots(e)[root_index]

"""
    reachable_variables(e::PathEdge)

This returns the internal data structure so **DO NOT MODIFY** the returned value unless you know what you are doing! Modification will corrupt the derivative graph. If you need to modify the result then copy first."""
reachable_variables(e::PathEdge) = e.reachable_variables

is_variable_reachable(e::PathEdge, variable_index::Integer) = reachable_variables(e)[variable_index]

"""
    to_tuple(e::PathEdge)

used for printing out a more readable version of Edge"""
to_tuple(e::PathEdge) = (e.top_vertex, e.bott_vertex)


num_reachable_variables(e::PathEdge) = sum(e.reachable_variables)

num_reachable_roots(e::PathEdge) = sum(e.reachable_roots)

num_uses(e::PathEdge) = num_reachable_variables(e) * num_reachable_roots(e)


can_delete(e::PathEdge) = !any(reachable_roots(e)) || !any(reachable_variables(e))

function Base.show(io::IO, a::PathEdge)
    print(io, "($(top_vertex(a)) $(bott_vertex(a))  $(num_uses(a)) $(value(a)) $(reachable_roots(a)) $(reachable_variables(a)))")
end

mask_roots!(e::PathEdge, mask::BitVector) = e.reachable_roots .= e.reachable_roots .& mask
mask_variables!(e::PathEdge, mask::BitVector) = e.reachable_variables .= e.reachable_variables .& mask

"""
    zero_roots!(e::PathEdge, root_mask::BitVector)

`root_mask` is a `BitVector` containing a 1 at index `i` if the edge is reachable from root `i` and a 0 otherwise. Changes the `reachable_roots` field of the edge to be 0 where `root_mask` is 1 and unchanged otherwise."""
zero_roots!(e::PathEdge, root_mask::BitVector) = @. e.reachable_roots = e.reachable_roots & (!root_mask)

"""
    zero_variables!(e::PathEdge, variable_mask::BitVector)

`variable_mask` is a `BitVector` containing a 1 at index `i` if the edge is reachable from variable `i` and a 0 otherwise. Changes the `reachable_variables` field of the edge to be 0 where `variable_mask` is 1 and unchanged otherwise."""
zero_variables!(e::PathEdge, variable_mask::BitVector) = @. e.reachable_variables = e.reachable_variables & (!variable_mask)


roots_resettable(e::PathEdge, root_mask::BitVector) = overlap(e.reachable_roots, root_mask)

variables_resettable(e::PathEdge, variable_mask::BitVector) = overlap(e.reachable_variables, variable_mask)

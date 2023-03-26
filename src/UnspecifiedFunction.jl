struct UnspecifiedFunction{V,D}
    name::Symbol
    variables::SVector{V,Node{SymbolicUtils.BasicSymbolic{Real},0}}
    derivatives::SVector{D,Node{SymbolicUtils.BasicSymbolic{Real},0}}

    UnspecifiedFunction(name::Symbol, variables::SVector{V,S}, derivatives::SVector{D,S}) where {D,V,S<:Node{SymbolicUtils.BasicSymbolic{Real},0}} = new{V,D}(name, variables, derivatives)
end
export UnspecifiedFunction

#special case for unspecified functions. Unfortunately this must go here rather than in ExpressionGraph.jl because this depends on Node which is defined in ExpressionGraph.jl
function derivative(uf::UnspecifiedFunction{V,D}, args::NTuple{N,Node{SymbolicUtils.BasicSymbolic{Real},0}}, ::Val{I}) where {I,V,D,N}
    new_derivs = MVector{D + 1,Node{SymbolicUtils.BasicSymbolic{Real},0}}(uf.derivatives..., args[I])
    return Node(UnspecifiedFunction(uf.name, uf.variables, new_derivs), MVector(uf.variables))
end


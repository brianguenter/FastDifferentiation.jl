
#TODO: probably want to add boolean operations so can sort Nodes.
# binary ops that return Bool
# for (f, Domain) in [(==) => Number, (!=) => Number,
#     (<=) => Real,   (>=) => Real,
#     (isless) => Real,
#     (<) => Real,   (> ) => Real,
#     (& ) => Bool,   (| ) => Bool,
#     xor => Bool]
# @eval begin
# promote_symtype(::$(typeof(f)), ::Type{<:$Domain}, ::Type{<:$Domain}) = Bool
# (::$(typeof(f)))(a::Symbolic{<:$Domain}, b::$Domain) = term($f, a, b, type=Bool)
# (::$(typeof(f)))(a::Symbolic{<:$Domain}, b::Symbolic{<:$Domain}) = term($f, a, b, type=Bool)
# (::$(typeof(f)))(a::$Domain, b::Symbolic{<:$Domain}) = term($f, a, b, type=Bool)
# end
# end

# struct Differential
#     expression::FastDifferentiation.Node #expression to take derivative of
#     with_respect_to::FastDifferentiation.Node #variable or internal node to take derivative wrt
# end

# derivative(f, args, v) = NoDeriv()

Base.:^(a::FastDifferentiation.Node, b::Integer) = simplify_check_cache(^, a, b, EXPRESSION_CACHE)

rules = Any[]

# Base.push!(a::Vector{T}, b::Number) where {T<:Node} = push!(a, Node(b)) #there should be a better way to do this.

Base.convert(::Type{Node}, a::T) where {T<:Real} = Node(a)
Base.promote_rule(::Type{<:Real}, ::Type{<:Node}) = Node

#convert two nodes with different number types to the correct new type. Various math functions will crash if this isn't defined.
function Base.promote_rule(a::Node{T,0}, b::Node{S,0}) where {T<:Real,S<:Real}
    ptype = promote_type(T, S)
    return Node(ptype(value(a))), Node(ptype(value(b)))
end

Base.conj(a::Node) = a #need to define this because dot and probably other linear algebra functions call this.
Base.adjoint(a::Node) = a
Base.transpose(a::Node) = a

# Pre-defined derivatives
import DiffRules
for (modu, fun, arity) ∈ DiffRules.diffrules(; filter_modules=(:Base, :SpecialFunctions, :NaNMath))
    fun in [:*, :+, :abs, :mod, :rem, :max, :min] && continue # special
    for i ∈ 1:arity

        expr = if arity == 1
            DiffRules.diffrule(modu, fun, :(args[1]))
        else
            DiffRules.diffrule(modu, fun, ntuple(k -> :(args[$k]), arity)...)[i]
        end

        push!(rules, expr)
        @eval derivative(::typeof($modu.$fun), args::NTuple{$arity,Any}, ::Val{$i}) = $expr
    end
end

#need to define because derivative functions can return inv
Base.inv(a::Node{typeof(/),2}) = children(a)[2] / children(a)[1]
Base.inv(a::Node) = 1 / a

#efficient explicit methods for most common cases
derivative(a::Node{T,1}, index::Val{1}) where {T} = derivative(value(a), (children(a)[1],), index)
derivative(a::Node{T,2}, index::Val{1}) where {T} = derivative(value(a), (children(a)[1], children(a)[2]), index)
derivative(a::Node{T,2}, index::Val{2}) where {T} = derivative(value(a), (children(a)[1], children(a)[2]), index)
derivative(a::Node, index::Val{i}) where {i} = derivative(value(a), (children(a)...,), index)
export derivative


function derivative(::typeof(*), args::NTuple{N,Any}, ::Val{I}) where {N,I}
    if N == 2
        return I == 1 ? args[2] : args[1]
    else
        return Node(*, deleteat!(collect(args), I)...) #simplify_check_cache will only be called for 2 arguments or less. Need to extedn to nary *, n> 2, if this is necessary.
    end
end

derivative(::typeof(+), args::NTuple{N,Any}, ::Val{I}) where {I,N} = Node(1)

# Special cases for leaf nodes with no children.
function derivative(a::Node{T,0}) where {T}
    if is_variable(a)
        return Node(1)
    else
        return Node(0)
    end
end



"""returns the leaf variables in a DAG. If a leaf is a Sym the assumption is that it is a variable. Leaves can also be numbers, which are not variables. Not certain how robust this is."""
variables(node::Node) = filter((x) -> is_variable(x), graph_leaves(node)) #SymbolicUtils changed, used to use SymbolicUtils.Sym for this test.


children(a::Node) = a.children


function Base.show(io::IO, a::Node)
    print(io, to_string(a))
end

function to_string(a::Node)
    function node_id(b::Node)
        # return "Node:$(b.node_value) id:$(objectid(b))"
        return "$(b.node_value)"
    end

    if arity(a) == 0
        return "$(node_id(a))"
    else
        if arity(a) == 1
            return "$(node_id(a))($(to_string(a.children[1])))"
        elseif arity(a) == 2
            return "($(to_string(a.children[1])) $(node_id(a)) $(to_string(a.children[2])))"
        else #Symbolics has expressions like +,* that can have any number of arguments, which translates to any number of Node children.
            return "($(a.node_value) $(foldl((x,y) -> x * " " * y,to_string.(a.children), init = "")))" #this is probably incredibly inefficient since it's O(n^2) in the length of the expression. Presumably there won't be many expresssions x+y+z+....., that are incredibly long. Best bet don't print out giant expressions.
        end
    end
end

function node_symbol(a::Node, variable_to_index::IdDict{Node,Int64})
    if is_tree(a)
        result = gensym() #create a symbol to represent the node
    elseif is_variable(a)
        result = :(input_variables[$(variable_to_index[a])])
    else
        result = value(a) #not a tree not a variable so is some kind of constant.
    end
    return result
end

"""Used to postorder function with multiple outputs"""
function postorder(roots::AbstractVector{T}) where {T<:Node}
    node_to_index = IdDict{Node,Int64}()
    nodes = Vector{Node}(undef, 0)
    variables = Vector{Node}(undef, 0)

    for root in roots
        _postorder_nodes!(root, nodes, variables, node_to_index)
    end
    return node_to_index, nodes, variables
end

"""returns vector of `Node` entries in the tree in postorder, i.e., if `result[i] == a::Node` then the postorder number of `a` is`i`. Not Multithread safe."""
function _postorder_nodes!(a::Node{T,N}, nodes::AbstractVector{S}, variables::AbstractVector{S}, visited::IdDict{Node,Int64}) where {T,N,S<:Node}
    if get(visited, a, nothing) === nothing
        if a.children !== nothing
            for child in a.children
                _postorder_nodes!(child, nodes, variables, visited)
            end
        elseif is_variable(a)
            push!(variables, a)
        end
        push!(nodes, a)
        visited[a] = length(keys(visited)) + 1 #node has not been added to visited yet so the count will be one less than needed.
    end
    return nothing
end

"""finds all the nodes in the graph and the number of times each node is visited in DFS."""
function all_nodes(a::Node, index_type=DefaultNodeIndexType)
    visited = IdDict{Node,index_type}()
    nodes = Vector{Node}(undef, 0)

    _all_nodes!(a, visited, nodes)
    return nodes
end

function _all_nodes!(node::Node, visited::IdDict{Node,T}, nodes::Vector{Node}) where {T<:Integer}
    tmp = get(visited, node, nothing)
    if tmp === nothing
        push!(nodes, node) #only add node to nodes once.
        if node.children !== nothing
            _all_nodes!.(node.children, Ref(visited), Ref(nodes))
        end
        visited[node] = 1
    else #already visited this node so don't have to recurse to children
        visited[node] += 1
    end

    return nothing
end

"""computes leaves of the graph `node`. This is inefficient since it allocates space to store all nodes and then searches through that vector to find the leaves."""
function graph_leaves(node::Node)
    result = Vector{Node}(undef, 0)
    nodes = all_nodes(node)

    for n in nodes
        if arity(n) === 0 #all nodes with no children are leaves.
            push!(result, n)
        end
    end

    return result
end

"""Returns an Array of variables with names corresponding to their indices in the Array. Example:
```julia
julia> make_variables(:x,3)
3-element Vector{FastDifferentiation.Node}:
 x1
 x2
 x3

julia> make_variables(:x,2,3)
2×3 Matrix{FastDifferentiation.Node}:
 x1_1  x1_2  x1_3
 x2_1  x2_2  x2_3

julia> make_variables(:x,2,3,2)
2×3×2 Array{FastDifferentiation.Node, 3}:
[:, :, 1] =
 x1_1_1  x1_2_1  x1_3_1
 x2_1_1  x2_2_1  x2_3_1

[:, :, 2] =
 x1_1_2  x1_2_2  x1_3_2
 x2_1_2  x2_2_2  x2_3_2
 ```
 """
function make_variables(name::Symbol, array_size::T...) where {T}
    result = Array{Node,length(array_size)}(undef, array_size...)

    for i in CartesianIndices((UnitRange.(1, array_size)))
        result[i] = Node(Symbol(name, join(i.I, "_")))
    end
    return result
end
export make_variables

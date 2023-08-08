

#until I can think of a better way of structuring the caching operation it will be a single global expression cache. This precludes multithreading, unfortunately. Many other parts of the algorithm are difficult to multithread. Processing is also so quick that only large graphs would benefit from multithreading. Don't know how common these will be.
const EXPRESSION_CACHE = IdDict()

function check_cache(a::Tuple{Vararg}, cache)
    cache_val = get(cache, a, nothing)
    if cache_val === nothing
        cache[a] = Node(a[1], a[2:end]...) #this should wrap everything, including basic numbers, in a Node object
    end

    return cache[a]
end

"""Clears the global expression cache. To maximize efficiency of expressions the differentation system automatically eliminates common subexpressions by checking for their existence in the global expression cache. Over time this cache can become arbitrarily large. Best practice is to clear the cache before you start defining expressions, define your expressions and then clear the cache."""
clear_cache() = empty!(EXPRESSION_CACHE)
export clear_cache

"""
@variables args...

create **FD** variables to use in symbolic expressions. Example:

```julia
julia> @variables x y
y

julia> f = x*y
(x * y)
```
"""
macro variables(args...)
    tmp = Expr(:block)
    for x in args
        push!(tmp.args, :($(esc(x)) = Node($(Meta.quot(x)))))
    end
    return tmp
end
export @variables

# #also add this inner constructor
#
#
# #end of code block to add

struct Node <: Number
    node_value
    children::Union{Nothing,MVector{N,Node}} where {N} #initially used SVector but this was incredibly inefficient. Possibly because the compiler was inlining the entire graph into a single static structure, which could lead to very long == and hashing times.


    Node(f::S, a) where {S} = new(f, MVector{1,Node}(Node(a)))
    Node(f::S, a, b) where {S} = new(f, MVector{2,Node}(Node(a), Node(b))) #if a,b not a Node convert them.

    Node(a::T) where {T<:Real} = new(a, nothing) #convert numbers to Node
    Node(a::T) where {T<:Node} = a #if a is already a special node leave it alone

    function Node(operation, args::MVector{N,T}) where {T<:Node,N} #use MVector rather than Vector. 40x faster.
        return new(operation, args)
    end

    Node(a::S) where {S<:Symbol} = new(a, nothing)
end



#convenience function to extract the fields from Node object to check cache
function check_cache(a::Node, cache)
    if children(a) !== nothing
        check_cache((value(a), children(a)...), cache)
    else
        check_cache((a,), cache)
    end
end


Base.zero(::Type{Node}) = Node(0)
Base.zero(::Node) = Node(0)
Base.one(::Type{Node}) = Node(1)
Base.one(::Node) = Node(1)

Broadcast.broadcastable(a::Node) = (a,)

value(a::Node) = a.node_value

arity(a::Node) = a.children === nothing ? 0 : length(a.children)


is_leaf(a::Node) = arity(a) == 0

is_tree(a::Node) = arity(a) >= 1


is_variable(a::Node) = isa(value(a), Symbol)


is_constant(a::Node) = !is_variable(a) && !is_tree(a)

is_negate(a::Node) = typeof(value(a)) == typeof(-) && arity(a) == 1
export is_negate

function constant_value(a::Node)
    if is_constant(a)
        return value(a)
    else
        return nothing
    end
end

error_message() = throw(ErrorException("FastDifferentiation.jl does not currently support comparison operations on FastDifferentiation expressions. Your code, or libraries called by your code, had a statement with a comparison operator such as x<y where either or both of x,y were FastDifferentiation expressions. You have to remove all these comparison operators for FastDifferentiation to work."))
Base.isless(::Node, ::Number) = error_message()
Base.isless(::Number, ::Node) = error_message()
Base.isless(::Node, ::Node) = error_message()

Base.iszero(a::Node) = value(a) == 0 #need this because sparse matrix and other code in linear algebra may call it. If it is not defined get a type promotion error.

function is_zero(a::Node)
    #this: value(a) == 0 would work but when add conditionals to the language if a is not a constant this will generate an expression graph instead of returning a bool value.
    if is_tree(a) || is_variable(a)
        return false
    elseif value(a) == 0
        return true
    else
        return false
    end
end


function is_one(a::Node)
    #this: value(a) == 1 would work but when add conditionals to the language if a is not a constant this will generate an expression graph instead of returning a bool value.
    if is_tree(a) || is_variable(a)
        return false
    elseif value(a) == 1
        return true
    else
        return false
    end
end


#Simple algebraic simplification rules for *,+,-,/. These are mostly safe, i.e., they will return exactly the same results as IEEE arithmetic. However multiplication by 0 always simplifies to 0, which is not true for IEEE arithmetic: 0*NaN=NaN, 0*Inf = NaN, for example. This should be a good tradeoff, since zeros are common in derivative expressions and can result in considerable expression simplification. Maybe later make this opt-out.

simplify_check_cache(a, b, c, cache) = check_cache((a, b, c), cache)

is_nary(a::Node) = arity(a) > 2
is_times(a::Node) = value(a) == *

is_nary_times(a::Node) = is_nary(a) && value(a) == typeof(*)



function simplify_check_cache(::typeof(^), a, b, cache)
    na = Node(a)
    nb = Node(b)
    if constant_value(na) !== nothing && constant_value(nb) !== nothing
        return Node(constant_value(na)^constant_value(nb))
    elseif value(nb) == 0 #IEEE-754 standard says any number, including Inf, NaN, etc., to the zero power is 1
        return Node(one(value(b))) #try to preserve number type
    elseif value(nb) == 1
        return a
    else
        return check_cache((^, na, nb), cache)
    end
end

function simplify_check_cache(::typeof(*), na, nb, cache)
    a = Node(na)
    b = Node(nb)

    #TODO sort variables so if y < x then x*y => y*x. The will automatically get commutativity.
    #c1*c2 = c3, (c1*x)*(c2*x) = c3*x
    if is_zero(a) && is_zero(b)
        return Node(value(a) + value(b)) #user may have mixed types for numbers so use automatic promotion to widen the type.
    elseif is_zero(a) #b is not zero
        return a #use this node rather than creating a zero since a has the type encoded in it
    elseif is_zero(b) #a is not zero
        return b #use this node rather than creating a zero since b has the type encoded in it
    elseif is_one(a)
        return b #At this point in processing the type of b may be impossible to determine, for example if b = sin(x) and the value of x won't be known till the expression is evaluated. No easy way to promote the type of b here if a has a wider type than b will eventually be determined to have. Example: a = BigFloat(1.0), b = sin(x). If the value of x is Float32 when the function is evaluated then would expect the type of the result to be BigFloat. But it will be Float32. Need to figure out a type of Node that will eventually generate code something like this: b = promote_type(a,b)(b) where the types of a,b will be known because this will be called in the generated Julia function for the derivative.
    elseif is_one(b)
        return a
    elseif is_constant(a) && value(a) == -1
        return -b
    elseif is_constant(b) && value(b) == -1
        return -a
    elseif constant_value(a) !== nothing && constant_value(b) !== nothing
        return Node(constant_value(a) * constant_value(b)) #this relies on the fact that objectid(Node(c)) == objectid(Node(c)) where c is a constant so don't have to check cache.
    elseif is_constant(b) #Move constant on right branch to left branch to simplify other simplification rules. Know that only one of a,b can be constant at this point so won't recurse. 
        return b * a
    elseif is_constant(a) && typeof(*) == typeof(value(b)) && is_constant(children(b)[1])
        return Node(value(children(b)[1]) * value(a)) * children(b)[2]
    elseif typeof(*) == typeof(value(a)) && typeof(*) == typeof(value(b)) && is_constant(children(b)[1]) && is_constant(children(a)[1])
        return Node(value(children(a)[1]) * value(children(b)[1])) * (children(b)[2] * children(a)[2])
    else
        return check_cache((*, a, b), cache)
    end
end

function simplify_check_cache(::typeof(+), na, nb, cache)
    a = Node(na)
    b = Node(nb)

    #TODO sort variables so if y < x then x*y => y*x. The will automatically get commutativity.
    #TODO add another check that moves all contants to the left and then does constant propagation
    #c1 + c2 = c3, (c1 + x)+(c2 + x) = c3+2*x

    if is_zero(a)
        return b
    elseif is_zero(b)
        return a
    elseif is_constant(a) && is_constant(b)
        return Node(value(a) + value(b))
    elseif is_constant(a) && typeof(+) == typeof(value(b)) && is_constant(children(b)[1]) #C1 + (C2 + x) => ((C1+C2)+x)
        return Node(value(a) + value(children(b)[1])) + children(b)[2]
    else
        return check_cache((+, a, b), cache)
    end
end

function simplify_check_cache(::typeof(/), na, nb, cache)
    a = Node(na)
    b = Node(nb)

    if is_one(b)
        return return a
    elseif is_constant(a) && is_constant(b)
        return Node(value(a) / value(b))
    else
        return check_cache((/, a, b), cache)
    end
end

function simplify_check_cache(::typeof(-), na, nb, cache)
    a = Node(na)
    b = Node(nb)
    if is_zero(b)
        return a
    elseif is_zero(a)
        return -b
    elseif is_constant(a) && is_constant(b)
        return Node(value(a) - value(b))
    else
        return check_cache((-, a, b), cache)
    end
end

simplify_check_cache(f::Any, na, cache) = check_cache((f, na), cache)

"""Special case only for unary -. No simplifications are currently applied to any other unary functions"""
function simplify_check_cache(::typeof(-), a, cache)
    na = Node(a) #this is safe because Node constructor is idempotent
    if arity(na) == 1 && typeof(value(na)) == typeof(-)
        return children(na)[1]
    elseif constant_value(na) !== nothing
        return Node(-value(na))
    elseif typeof(*) == typeof(value(na)) && constant_value(children(na)[1]) !== nothing
        return Node(-value(children(na)[1])) * children(na)[2]
    else
        return check_cache((-, na), cache)
    end
end

Base.:^(a::FastDifferentiation.Node, b::Integer) = simplify_check_cache(^, a, b, EXPRESSION_CACHE)

rules = Any[]

Base.convert(::Type{Node}, a::T) where {T<:Real} = Node(a)
Base.promote_rule(::Type{<:Real}, ::Type{Node}) = Node

function Base.:-(a::AbstractArray{<:Node,N}) where {N}
    if length(a) == 0
        return a
    else
        return .-(a)
    end
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

function Base.inv(a::Node)
    if typeof(value(a)) === /
        return children(a)[2] / children(a)[1]
    else
        return 1 / a
    end
end


#need special case for sincos because it returns a 2 tuple. Also Diffrules.jl does not define a differentiation rule for sincos.
Base.sincos(x::Node) = (sin(x), cos(x)) #this will be less efficient than sincos. TODO figure out a better way.

derivative(::typeof(abs), arg::Tuple{T}, ::Val{1}) where {T} = arg[1] / abs(arg[1])

function derivative(::typeof(*), args::NTuple{N,Any}, ::Val{I}) where {N,I}
    if N == 2
        return I == 1 ? args[2] : args[1]
    else
        return Node(*, deleteat!(collect(args), I)...) #TODO: simplify_check_cache will only be called for 2 arguments or less. Need to extedn to nary *, n> 2, if this is necessary.
    end
end

derivative(::typeof(+), args::NTuple{N,Any}, ::Val{I}) where {I,N} = Node(1)

function derivative(a::Node, index::Val{1})
    if arity(a) == 1
        return derivative(value(a), (children(a)[1],), index)
    elseif arity(a) == 2
        return derivative(value(a), (children(a)[1], children(a)[2]), index)
    else
        throw(ErrorException("should never get here"))
    end
end

function derivative(a::Node, index::Val{2})
    if arity(a) == 2
        return derivative(value(a), (children(a)[1], children(a)[2]), index)
    else
        throw(ErrorException("should never get here"))
    end
end

derivative(a::Node, index::Val{i}) where {i} = derivative(value(a), (children(a)...,), index)

# Special cases for leaf nodes with no children.
function derivative(a::Node)
    if is_leaf(a)
        if is_variable(a)
            return one(Node)
        else
            return zero(Node)
        end
    else
        throw(ErrorException("should never reach this point in the code"))
    end
end
export derivative



"""returns the leaf variables in a DAG. If a leaf is a Sym the assumption is that it is a variable. Leaves can also be numbers, which are not variables. Not certain how robust this is."""
variables(node::Node) = filter((x) -> is_variable(x), graph_leaves(node))

function variables(a::AbstractArray{T}) where {T<:Node}
    visited = IdDict{Node,Int64}()
    for root in a
        all_nodes!(root, visited)
    end
    return filter((x) -> is_variable(x), collect(keys(visited)))
end

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
function _postorder_nodes!(a::Node, nodes::AbstractVector{S}, variables::AbstractVector{S}, visited::IdDict{Node,Int64}) where {S<:Node}
    if get(visited, a, nothing) === nothing
        if arity(a) != 0
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

    all_nodes!(a, visited)
    return collect(keys(visited))
    return nodes
end

function all_nodes!(node::N, visited::IdDict{Node,T}) where {T<:Integer,N<:Node}
    tmp = get(visited, node, nothing)
    if tmp === nothing
        visited[node] = 1
        if arity(node) != 0
            all_nodes!.(node.children, Ref(visited))
        end
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

"""
make_variables(name::Symbol,array_size::T...)

Returns an Array of variables with names corresponding to their indices in the Array. Example:
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


#create methods that accept Node arguments for all mathematical functions.
@number_methods(Node, simplify_check_cache(f, a, EXPRESSION_CACHE), simplify_check_cache(f, a, b, EXPRESSION_CACHE)) #create methods for standard functions that take Node instead of Number arguments. Check cache to see if these arguments have been seen before.

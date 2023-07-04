

#until I can think of a better way of structuring the caching operation it will be a single global expression cache. This precludes multithreading, unfortunately. Many other parts of the algorithm are difficult to multithread. Processing is also so quick that only large graphs would benefit from multithreading. Don't know how common these will be.
EXPRESSION_CACHE = IdDict()

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

struct Node{T,N} <: Number
    node_value::T
    children::Union{MVector{N,Node},Nothing} #initially used SVector but this was incredibly inefficient. Possibly because the compiler was inlining the entire graph into a single static structure, which could lead to very long == and hashing times.

    Node(f::S, a) where {S} = new{S,1}(f, MVector{1}(Node(a)))
    Node(f::S, a, b) where {S} = new{S,2}(f, MVector{2}(Node(a), Node(b))) #if a,b not a Node convert them.

    Node(a::T) where {T<:Real} = new{T,0}(a, nothing) #convert numbers to Node
    Node(a::T) where {T<:Node} = a #if a is already a special node leave it alone

    Node(a::AutomaticDifferentiation.NoDeriv) = new{AutomaticDifferentiation.NoDeriv,0}(a, nothing) #TODO: this doesn't seem like it should ever be called.

    function Node(operation, args::MVector{N,T}) where {T<:Node,N} #use MVector rather than Vector. 40x faster.
        ntype = typeof(operation)
        return new{ntype,N}(operation, args)
    end

    Node(a::S) where {S<:Symbol} = new{S,0}(a, nothing)
end

#convenience function to extract the fields from Node object to check cache
function check_cache(a::Node{T,N}, cache) where {T,N}
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

arity(::Node{T,N}) where {T,N} = N


is_leaf(::Node{T,0}) where {T} = true
is_leaf(::Node{T,N}) where {T,N} = false

is_tree(::Node{T,N}) where {T,N} = N >= 1


is_variable(a::Node) = isa(value(a), Symbol)


is_constant(a::Node) = !is_variable(a) && !is_tree(a)


function constant_value(a::Node)
    if is_constant(a)
        return value(a)
    else
        return nothing
    end
end

Base.iszero(a::Node) = value(a) == 0 #need this because sparse matrix and other code in linear algebra may call it. If it is not defined get a type promotion error.
function is_zero(a::Node)
    if is_tree(a) || is_variable(a)
        return false
    elseif value(a) == 0
        return true
    else
        return false
    end
end


function is_one(a::Node)
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

is_nary(a::Node{T,N}) where {T,N} = N > 2
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
        return Node(value(children(a)[1] * value(children(b)[1]))) * (children(b)[2] * children(a)[2])
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
    else
        return check_cache((-, na), cache)
    end
end

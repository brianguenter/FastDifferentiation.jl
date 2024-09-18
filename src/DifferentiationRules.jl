# Pre-defined derivatives
import DiffRules
#Special case rules for rules in DiffRules that use ?: or if...else, neither of which will work when add conditionals

#airybix and airyprimex don't work in FastDifferentiation. airybix(x) where x is a Node causes a stack overflow. So no diffrule defined, although airybix uses if...else
#LogExpFunctions.xlogy uses ?: so doesn't work with FastDifferentiation. Most of the functions in this package don't work with FastDifferentiation.

DiffRules.@define_diffrule Base.:^(x, y) = :($y * ($x^($y - 1))), :(if_else($x isa Real && $x <= 0, Base.oftype(float($x), NaN), ($x^$y) * log($x)))


for (modu, fun, arity) ∈ DiffRules.diffrules(; filter_modules=(:Base, :SpecialFunctions, :NaNMath))
    fun in [:*, :+, :abs, :mod, :rem, :max, :min] && continue # special
    for i ∈ 1:arity

        expr = if arity == 1
            DiffRules.diffrule(modu, fun, :(args[1]))
        else
            DiffRules.diffrule(modu, fun, ntuple(k -> :(args[$k]), arity)...)[i]
        end

        @eval derivative(::typeof($modu.$fun), args::NTuple{$arity,Any}, ::Val{$i}) = $expr
    end
end

derivative(::typeof(abs), arg::Tuple{T}, ::Val{1}) where {T} = arg[1] / abs(arg[1])

function derivative(::typeof(*), args::NTuple{N,Any}, ::Val{I}) where {N,I}
    if N == 2
        return I == 1 ? args[2] : args[1]
    else
        return Node(*, deleteat!(collect(args), I)...) #TODO: simplify_check_cache will only be called for 2 arguments or less. Need to extend to nary *, n> 2, if this is necessary.
    end
end

derivative(::typeof(+), args::NTuple{N,Any}, ::Val{I}) where {I,N} = Node(1)
derivative(::NoOp, arg::Tuple{T}, ::Val{1}) where {T} = 1.0


function_variable_derivative(a::Node, index::Val{i}) where {i} = check_cache((Differential, children(a)[i]))

# These functions are primarily used to do error checking on expressions
function derivative(a::Node, index::Val{1})
    if is_conditional(a)
        throw(conditional_error(a))
    elseif is_variable_function(a)
        return function_variable_derivative(a, index)
    elseif arity(a) == 1
        return derivative(value(a), (children(a)[1],), index)
    elseif arity(a) == 2
        return derivative(value(a), (children(a)[1], children(a)[2]), index)
    else
        throw(ErrorException("should never get here"))
    end
end

function derivative(a::Node, index::Val{2})
    if is_conditional(a)
        throw(conditional_error(a))
    elseif is_variable_function(a)
        return function_variable_derivative(a, index)
    elseif arity(a) == 2
        return derivative(value(a), (children(a)[1], children(a)[2]), index)
    else
        throw(ErrorException("should never get here"))
    end
end

function derivative(a::Node, index::Val{i}) where {i}
    if is_conditional(a)
        throw(conditional_error(a))
    elseif is_variable_function(a)
        return function_variable_derivative(a, index)
    else
        return derivative(value(a), (children(a)...,), index)
    end
end
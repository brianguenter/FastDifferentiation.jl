abstract type SymbolicWrapper end
abstract type RuntimeWrapper end
struct NoConditionals <: SymbolicWrapper
    cached_derivative::DerivativeGraph
end

struct Conditionals <: SymbolicWrapper
    cached_derivatives::Dict{Vector{Bool},DerivativeGraph}
    root_functions::Vector{Bool}
    conditional_function::RuntimeGeneratedFunction

    function Conditionals(roots, input_variables)
        cond_func = make_conditional_function(roots, input_variables)
        cached_deriv = Dict{Vector{Bool},DerivativeGraph}(undef, 0)
        return new(cached_deriv, roots, cond_func)
    end
end

struct NoConditionalsRuntime <: RuntimeWrapper
    cached_function::RuntimeGeneratedFunction
end

struct ConditionalsRuntime <: RuntimeWrapper
    dgraph_cache::Conditionals
    cached_functions::Dict{DerivativeGraph,RuntimeGeneratedFunction}
end

#make_function returns these instead of a RuntimeGeneratedFunction
(a::NoConditionalsRuntime)(input_values::T) where {T<:Real} = a.cached_derivative(input_values)
function (a::ConditionalsRuntime)(input_values::T) where {T<:Real}
    tmp = a.dgraph_cache
    a.cached_derivatives[a.conditional_function(input_values)](input_values)

"""Creates a cache object for derivative graphs with conditionals. To evaluate the derivative for inputs `𝐱` you use the object as a function. Example:

```@variables x y

    gr = dgraph([x^2,y^3],[x,y])
    jacobian(gr)"""
function dgraph(roots::AbstractVector{Node}, input_variables::AbstractVector{Node})
    postorder_number = IdDict{Node,index_type}()

    (postorder_number, nodes, _) = postorder(roots)

    if length(ifelse_nodes(nodes)) == 0 && length(cond_nodes(nodes)) == 0 #if no conditional boolean values passed in then the graph must not contain boolean or ifelse nodes
        return NoConditionals(DerivativeGraph(roots))
    else
        return Conditionals(roots, input_variables)
    end
end


derivative_graph(wrapper::NoConditionals) = wrapper.cached_derivative
function derivative_graph!(wrapper::Conditionals, input_values::T) where {T<:Real}
    conds = wrapper.conditional_function(input_values)
    tmp = get_index(wrapper.cached_derivatives, conds, nothing)

    if tmp !== nothing
        return tmp
    else
        gr = DerivativeGraph(instantiate_conditionals(wrapper.root_functions, conds))
        wrapper.cached_derivatives[conds] = gr
        return gr
    end
end

"""Traverse the dag defined by `func` """
function instantiate_conditionals(func::Vector{Node}, conditional_values::AbstractVector{T}) where {T<:Bool}
end

function make_conditional_function(roots, input_variables)
end

is_conditional(a::Node) = is_tree(a) && value(a) in BOOL_OPS

ifelse_nodes(nodes::Vector{Node}) = filter(x -> x === ifelse, nodes) #nodes(gr) is sorted by postorder number for results are as well.

bool_nodes(nodes::Vector{Node}) = filter(x -> x in BOOL_OPS, nodes) #nodes(gr) is sorted by postorder number for results are as well. 



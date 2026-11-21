# Prototype: a tracing JIT for FastDifferentiation derivatives.
#
# Idea (from LuaJIT / Dynamo / Dr.Jit): run the function on the actual input values, record the straight line
# code that was executed as an expression graph (FD's Node DAG is already an SSA form with common
# subexpression sharing), and log a *guard* for every branch that was taken. The recorded trace has no branches,
# so FD's ordinary `jacobian` applies. Each trace is compiled into ONE function that checks its guards (returning
# `nothing` if one fails) and otherwise computes the Jacobian. Traces are cached (most recently used first, with
# an LRU cap) and a new one is recorded only when no cached trace's guards hold.
#
# Nothing here reads Julia's IR: because `Tracked` carries a concrete value, `x < y` returns a real Bool and
# Julia's own `if`, `while`, `&&` and `throw` just work.

module TracingPrototype

using FastDifferentiation
using LinearAlgebra
using Printf

const FD = FastDifferentiation
const Node = FD.Node
const RGF = FD.RuntimeGeneratedFunctions # FD already depends on it

export Tracked, TracingJacobian, eval_jacobian, eval_jacobian!, loop_jacobian, summary_line, trace_report

# ------------------------------------------------------------------ recording context

const BIT_CAPACITY = 2048 # decisions. The buffers double in the rare case that a function makes more.
const INT_CAPACITY = 256

"""
Every comparison outcome (a bit, packed 64 to a word) and every integer taken from a traced value, in execution
order. The buffers are preallocated and written by index with a counter, not appended to; they only grow if a
function makes more decisions than fit.
"""
mutable struct DecisionLog
    bits::BitVector
    ints::Vector{Int}
    nbits::Int
    nints::Int
end
DecisionLog() = DecisionLog(falses(BIT_CAPACITY), Vector{Int}(undef, INT_CAPACITY), 0, 0)
reset!(l::DecisionLog) = (l.nbits = 0; l.nints = 0; l)

@noinline grow_bits!(l::DecisionLog) = (resize!(l.bits, 2 * length(l.bits)); nothing)
@noinline grow_ints!(l::DecisionLog) = (resize!(l.ints, 2 * length(l.ints)); nothing)

@inline function log_bit!(l::DecisionLog, outcome::Bool)
    n = l.nbits + 1
    n > length(l.bits) && grow_bits!(l)
    @inbounds l.bits[n] = outcome
    l.nbits = n
    return outcome
end
@inline function log_int!(l::DecisionLog, k::Integer)
    n = l.nints + 1
    n > length(l.ints) && grow_ints!(l)
    @inbounds l.ints[n] = Int(k)
    l.nints = n
    return k
end

@inline mix(h::UInt64, w::UInt64) = (h = (h ⊻ w) * 0x00000100000001b3; h ⊻ (h >> 32))

# A hash of the decisions made so far: whole 64 bit words of the packed bits, then the integers.
function signature(l::DecisionLog)
    h = 0xcbf29ce484222325 ⊻ UInt64(l.nbits) ⊻ (UInt64(l.nints) << 32)
    chunks = l.bits.chunks
    nfull = l.nbits >> 6
    @inbounds for i in 1:nfull
        h = mix(h, chunks[i])
    end
    r = l.nbits & 63
    if r != 0
        @inbounds h = mix(h, chunks[nfull+1] & ((UInt64(1) << r) - 1)) # the bits past the end of the log are stale
    end
    @inbounds for i in 1:l.nints
        h = mix(h, reinterpret(UInt64, l.ints[i]))
    end
    return h
end

mutable struct Recording
    guards::Vector{Pair{Node,Bool}} # condition expression => outcome that was observed (for guards compiled into the trace)
    log::DecisionLog                # every decision, in order (for replay)
end
Recording() = Recording(Pair{Node,Bool}[], DecisionLog())
const CURRENT = Ref{Union{Nothing,Recording}}(nothing)

function current_recording()
    rec = CURRENT[]
    rec === nothing && error("Tracked values can only be compared while a trace is being recorded")
    return rec
end

# A guard to be compiled into the trace: this condition must come out with this outcome.
function node_guard!(cond::Node, outcome::Bool)
    rec = current_recording()
    any(g -> g.first === cond, rec.guards) || push!(rec.guards, cond => outcome)
    return outcome
end

# A comparison or other Bool test: log the outcome (for replay) and the guard (for compiled guards).
# Constants are logged too, because the value-only replay number cannot tell a constant from a traced value, and the
# replay has to see exactly the same sequence of decisions.
function branch!(cond::Node, outcome::Bool, constant::Bool)
    rec = current_recording()
    log_bit!(rec.log, outcome)
    constant || node_guard!(cond, outcome)
    return outcome
end

# An integer taken from a traced value (index, loop bound, floor, ...).
function concrete!(k::Integer)
    log_int!(current_recording().log, k)
    return k
end

# ------------------------------------------------------------------ the traced number

"""A number that carries its actual value and the expression graph that computed it."""
struct Tracked <: Real
    val::Float64
    node::Node
end
Tracked(x::Tracked) = x
Tracked(x::Real) = Tracked(Float64(x), Node(x))

Base.promote_rule(::Type{Tracked}, ::Type{<:Real}) = Tracked
Base.convert(::Type{Tracked}, x::Real) = Tracked(x)
Base.convert(::Type{Tracked}, x::Tracked) = x
Base.zero(::Type{Tracked}) = Tracked(0.0)
Base.zero(::Tracked) = Tracked(0.0)
Base.one(::Type{Tracked}) = Tracked(1.0)
Base.one(::Tracked) = Tracked(1.0)
Base.float(x::Tracked) = x
Base.show(io::IO, x::Tracked) = print(io, "Tracked(", x.val, ")")

for f in (:sin, :cos, :tan, :exp, :log, :sqrt, :abs, :tanh, :asin, :acos, :atan, :sinh, :cosh)
    @eval Base.$f(a::Tracked) = Tracked(Base.$f(a.val), Base.$f(a.node))
end
Base.:-(a::Tracked) = Tracked(-a.val, -a.node)
Base.:+(a::Tracked) = a
for f in (:+, :-, :*, :/, :^, :atan, :hypot)
    @eval Base.$f(a::Tracked, b::Tracked) = Tracked(Base.$f(a.val, b.val), Base.$f(a.node, b.node))
    @eval Base.$f(a::Tracked, b::Real) = Base.$f(a, Tracked(b))
    @eval Base.$f(a::Real, b::Tracked) = Base.$f(Tracked(a), b)
end
Base.:^(a::Tracked, b::Integer) = Tracked(a.val^b, a.node^b)

# Comparisons return the real Bool and record a guard. This is what lets ordinary Julia control flow run.
for op in (:<, :>, :<=, :>=, :(==))
    @eval function Base.$op(a::Tracked, b::Tracked)
        r = Base.$op(a.val, b.val)
        if FD.is_constant(a.node) && FD.is_constant(b.node)
            return branch!(a.node, r, true)
        end
        return branch!(Base.$op(a.node, b.node), r, false)
    end
    @eval Base.$op(a::Tracked, b::Real) = Base.$op(a, Tracked(b))
    @eval Base.$op(a::Real, b::Tracked) = Base.$op(Tracked(a), b)
end
Base.isless(a::Tracked, b::Tracked) = a < b # so that Base's max, min, clamp, sort, ifelse work
Base.isless(a::Tracked, b::Real) = a < b
Base.isless(a::Real, b::Tracked) = a < b
Base.isless(a::Tracked, b::AbstractFloat) = a < b
Base.isless(a::AbstractFloat, b::Tracked) = a < b

# Using a traced value where Julia needs an integer (indices, loop bounds). The guard says the value must be the same.
function (::Type{T})(a::Tracked) where {T<:Integer}
    k = T(a.val) # InexactError if the value is not an integer
    FD.is_constant(a.node) || node_guard!(a.node == a.val, true)
    return concrete!(k)
end
Base.to_index(a::Tracked) = Int(a)

function Base.floor(::Type{T}, a::Tracked) where {T<:Integer}
    k = floor(T, a.val)
    if !FD.is_constant(a.node)
        node_guard!(a.node >= k, true)
        node_guard!(a.node < k + 1, true)
    end
    return concrete!(k)
end
Base.floor(a::Tracked) = Tracked(float(floor(Int, a)))

# mod is a node of its own, so a guard on Int(mod(u, 4)) is "mod(u, 4) == 1" and holds for u = 1, 5, 9, ...
Base.mod(a::Tracked, b::Real) = Tracked(mod(a.val, b), mod(a.node, b))
# FD on this branch has no derivative rules for mod (they were added on the TrueConditionals branch). The same two rules:
FD.derivative(::typeof(mod), args::Tuple{Any,Any}, ::Val{1}) = Node(1)
FD.derivative(::typeof(mod), args::Tuple{Any,Any}, ::Val{2}) = -(args[1] - mod(args[1], args[2])) / args[2]
Base.rem(a::Tracked, b::Real) = Tracked(rem(a.val, b), rem(a.node, b))
FD.derivative(::typeof(rem), args::Tuple{Any,Any}, ::Val{1}) = Node(1)
FD.derivative(::typeof(rem), args::Tuple{Any,Any}, ::Val{2}) = -(args[1] - rem(args[1], args[2])) / args[2]
Base.ldexp(a::Tracked, n::Integer) = a * 2.0^n

# ------------------------------------------------------------------ more of the numeric interface that packages use

Base.eps(::Type{Tracked}) = eps(Float64)
Base.eps(a::Tracked) = eps(a.val)
Base.typemax(::Type{Tracked}) = Tracked(Inf)
Base.typemin(::Type{Tracked}) = Tracked(-Inf)
Base.floatmin(::Type{Tracked}) = Tracked(floatmin(Float64))
Base.floatmax(::Type{Tracked}) = Tracked(floatmax(Float64))
Base.precision(::Type{Tracked}) = precision(Float64)
Base.precision(::Tracked) = precision(Float64)
Base.hash(a::Tracked, h::UInt) = hash(a.val, h)
Base.isequal(a::Tracked, b::Tracked) = a == b

# tests that return a Bool record a guard, like comparisons
for f in (:isnan, :isinf, :isfinite, :signbit)
    @eval function Base.$f(a::Tracked)
        r = Base.$f(a.val)
        FD.is_constant(a.node) && return branch!(a.node, r, true)
        return branch!(Base.$f(a.node), r, false)
    end
end
Base.sign(a::Tracked) = a < 0 ? Tracked(-1.0) : (a > 0 ? Tracked(1.0) : Tracked(0.0))
Base.copysign(a::Tracked, b::Tracked) = signbit(b) ? -abs(a) : abs(a)
Base.copysign(a::Tracked, b::Real) = signbit(b) ? -abs(a) : abs(a)
Base.copysign(a::Real, b::Tracked) = signbit(b) ? -abs(a) : abs(a)

for f in (:log2, :log10, :log1p, :expm1, :cbrt, :exp2, :exp10, :asinh, :acosh, :atanh, :sec, :csc, :cot)
    @eval Base.$f(a::Tracked) = Tracked(Base.$f(a.val), Base.$f(a.node))
end
Base.sincos(a::Tracked) = (sin(a), cos(a))

# rounding to an integer: the value is piecewise constant, so the guards say which interval the input is in
function Base.round(::Type{T}, a::Tracked) where {T<:Integer}
    k = round(T, a.val)
    if !FD.is_constant(a.node)
        node_guard!(a.node >= k - 0.5, true)
        node_guard!(a.node <= k + 0.5, true)
    end
    return concrete!(k)
end
function Base.ceil(::Type{T}, a::Tracked) where {T<:Integer}
    k = ceil(T, a.val)
    if !FD.is_constant(a.node)
        node_guard!(a.node > k - 1, true)
        node_guard!(a.node <= k, true)
    end
    return concrete!(k)
end
function Base.trunc(::Type{T}, a::Tracked) where {T<:Integer}
    return a.val >= 0 ? floor(T, a) : ceil(T, a)
end
Base.unsafe_trunc(::Type{T}, a::Tracked) where {T<:Integer} = trunc(T, a)
# nextfloat/prevfloat are used on constants such as thresholds
function Base.nextfloat(a::Tracked)
    FD.is_constant(a.node) || error("nextfloat of a Tracked value that is not a constant")
    return Tracked(nextfloat(a.val))
end
function Base.prevfloat(a::Tracked)
    FD.is_constant(a.node) || error("prevfloat of a Tracked value that is not a constant")
    return Tracked(prevfloat(a.val))
end
Base.rem2pi(a::Tracked, ::RoundingMode{:Nearest}) = a - 2π * round(Int, a / 2π)
function Base.exponent(a::Tracked)
    k = exponent(a.val)
    if !FD.is_constant(a.node)
        node_guard!(abs(a.node) >= 2.0^k, true)
        node_guard!(abs(a.node) < 2.0^(k + 1), true)
    end
    return concrete!(k)
end
Base.round(a::Tracked) = Tracked(float(round(Int, a)))
Base.ceil(a::Tracked) = Tracked(float(ceil(Int, a)))
Base.trunc(a::Tracked) = Tracked(float(trunc(Int, a)))

# ------------------------------------------------------------------ replay: find the trace for a new input by running f once

# Replaying runs f once on Probe numbers, which write their decisions into this log (by index, no allocation and no
# exceptions). The cached trace with exactly the same decisions is then looked up. A different path simply gives a
# decision sequence that no cached trace has.
const REPLAY_LOG = DecisionLog() # not thread safe, like FD's expression cache

@inline replay_bit(outcome::Bool) = log_bit!(REPLAY_LOG, outcome)
@inline replay_int(k::Integer) = log_int!(REPLAY_LOG, k)

"""A number that carries only its value. Every decision it takes is written to the replay log."""
struct Probe <: Real
    v::Float64
end
Probe(x::Probe) = x
Probe(x::Real) = Probe(Float64(x))

Base.promote_rule(::Type{Probe}, ::Type{<:Real}) = Probe
Base.convert(::Type{Probe}, x::Real) = Probe(x)
Base.convert(::Type{Probe}, x::Probe) = x
Base.zero(::Type{Probe}) = Probe(0.0)
Base.zero(::Probe) = Probe(0.0)
Base.one(::Type{Probe}) = Probe(1.0)
Base.one(::Probe) = Probe(1.0)
Base.float(x::Probe) = x
Base.show(io::IO, x::Probe) = print(io, "Probe(", x.v, ")")

for f in (:sin, :cos, :tan, :exp, :log, :sqrt, :abs, :tanh, :asin, :acos, :atan, :sinh, :cosh,
    :log2, :log10, :log1p, :expm1, :cbrt, :exp2, :exp10, :asinh, :acosh, :atanh, :sec, :csc, :cot)
    @eval Base.$f(a::Probe) = Probe(Base.$f(a.v))
end
Base.:-(a::Probe) = Probe(-a.v)
Base.:+(a::Probe) = a
for f in (:+, :-, :*, :/, :^, :atan, :hypot)
    @eval Base.$f(a::Probe, b::Probe) = Probe(Base.$f(a.v, b.v))
    @eval Base.$f(a::Probe, b::Real) = Base.$f(a, Probe(b))
    @eval Base.$f(a::Real, b::Probe) = Base.$f(Probe(a), b)
end
Base.:^(a::Probe, b::Integer) = Probe(a.v^b)
for op in (:<, :>, :<=, :>=, :(==))
    @eval Base.$op(a::Probe, b::Probe) = replay_bit(Base.$op(a.v, b.v))
    @eval Base.$op(a::Probe, b::Real) = Base.$op(a, Probe(b))
    @eval Base.$op(a::Real, b::Probe) = Base.$op(Probe(a), b)
end
Base.isless(a::Probe, b::Probe) = a < b
Base.isless(a::Probe, b::Real) = a < b
Base.isless(a::Real, b::Probe) = a < b
Base.isless(a::Probe, b::AbstractFloat) = a < b
Base.isless(a::AbstractFloat, b::Probe) = a < b
Base.isequal(a::Probe, b::Probe) = a == b
for f in (:isnan, :isinf, :isfinite, :signbit)
    @eval Base.$f(a::Probe) = replay_bit(Base.$f(a.v))
end
Base.sign(a::Probe) = a < 0 ? Probe(-1.0) : (a > 0 ? Probe(1.0) : Probe(0.0))
Base.copysign(a::Probe, b::Probe) = signbit(b) ? -abs(a) : abs(a)
Base.copysign(a::Probe, b::Real) = signbit(b) ? -abs(a) : abs(a)
Base.copysign(a::Real, b::Probe) = signbit(b) ? -abs(a) : abs(a)
Base.sincos(a::Probe) = (sin(a), cos(a))
Base.mod(a::Probe, b::Real) = Probe(mod(a.v, b))
Base.rem(a::Probe, b::Real) = Probe(rem(a.v, b))
Base.ldexp(a::Probe, n::Integer) = Probe(ldexp(a.v, n))

(::Type{T})(a::Probe) where {T<:Integer} = replay_int(T(a.v))
Base.to_index(a::Probe) = Int(a)
Base.floor(::Type{T}, a::Probe) where {T<:Integer} = replay_int(floor(T, a.v))
Base.round(::Type{T}, a::Probe) where {T<:Integer} = replay_int(round(T, a.v))
Base.ceil(::Type{T}, a::Probe) where {T<:Integer} = replay_int(ceil(T, a.v))
Base.trunc(::Type{T}, a::Probe) where {T<:Integer} = a.v >= 0 ? floor(T, a) : ceil(T, a)
Base.unsafe_trunc(::Type{T}, a::Probe) where {T<:Integer} = trunc(T, a)
Base.floor(a::Probe) = Probe(float(floor(Int, a)))
Base.round(a::Probe) = Probe(float(round(Int, a)))
Base.ceil(a::Probe) = Probe(float(ceil(Int, a)))
Base.trunc(a::Probe) = Probe(float(trunc(Int, a)))
Base.rem2pi(a::Probe, ::RoundingMode{:Nearest}) = a - 2π * round(Int, a / 2π)
Base.exponent(a::Probe) = replay_int(exponent(a.v))
Base.nextfloat(a::Probe) = Probe(nextfloat(a.v))
Base.prevfloat(a::Probe) = Probe(prevfloat(a.v))
Base.eps(::Type{Probe}) = eps(Float64)
Base.eps(a::Probe) = eps(a.v)
Base.typemax(::Type{Probe}) = Probe(Inf)
Base.typemin(::Type{Probe}) = Probe(-Inf)
Base.floatmin(::Type{Probe}) = Probe(floatmin(Float64))
Base.floatmax(::Type{Probe}) = Probe(floatmax(Float64))
Base.precision(::Type{Probe}) = precision(Float64)
Base.precision(::Probe) = precision(Float64)
Base.hash(a::Probe, h::UInt) = hash(a.v, h)

"""
Runs `f` on Probe numbers built from a flat input vector. The types of `f`, of the argument buffers and of the ranges are
all type parameters, so the argument tuple handed to `f` is concretely typed: no boxing, no `Vector{Any}`, no dynamic splat.
"""
struct Replayer{F,B<:Tuple,R<:Tuple}
    f::F
    bufs::B   # for each argument: `nothing` for a number, or a preallocated Array{Probe}
    ranges::R # where each argument sits in the flat input vector
end

function make_replayer(f, args)
    bufs = Tuple(a isa Number ? nothing : Array{Probe}(undef, size(a)) for a in args)
    return Replayer(f, bufs, Tuple(arg_ranges(args)))
end

@inline probe_arg(::Nothing, rng, x) = Probe(x[first(rng)])
@inline function probe_arg(buf::AbstractArray, rng, x)
    j = 0
    for idx in rng
        j += 1
        buf[j] = Probe(x[idx])
    end
    return buf
end

function run_replayer(rp::Replayer, x)
    args = map((b, r) -> probe_arg(b, r, x), rp.bufs, rp.ranges)
    return rp.f(args...)
end

# ------------------------------------------------------------------ flattening arguments

flatten(args) = Float64[v for a in args for v in (a isa Number ? (a,) : vec(a))]
function arg_ranges(args)
    ranges = UnitRange{Int}[]
    pos = 0
    for a in args
        len = a isa Number ? 1 : length(a)
        push!(ranges, (pos+1):(pos+len))
        pos += len
    end
    return ranges
end
function tracked_args(args, ranges, vals, vars)
    return Any[a isa Number ? Tracked(vals[first(r)], vars[first(r)]) :
               reshape([Tracked(vals[i], vars[i]) for i in r], size(a)) for (a, r) in zip(args, ranges)]
end

as_nodes(o::Tracked) = Node[o.node]
as_nodes(o::Real) = Node[Node(o)]
as_nodes(o::AbstractArray) = Node[x isa Tracked ? x.node : Node(x) for x in vec(o)]

# ------------------------------------------------------------------ one function per trace: guards, then the Jacobian

"""
Generate one function of the input vector that first checks the guards and then computes the Jacobian. Guard
conditions and Jacobian entries share subexpressions, so the forward computation is done once.

The calling convention is fixed when the function is generated:
  - `in_place=false`: `fn(x)` returns a newly allocated Jacobian, or `nothing` if a guard fails.
  - `in_place=true`:  `fn(result, x)` writes the Jacobian into `result` and returns `true`, or returns `false` if a
    guard fails. The guards run before anything is written, so a failed guard leaves `result` untouched and the
    caller can try another trace with the same array. `init_with_zeros=true` emits a fill of `result` with zeros
    (only if the Jacobian has structural zeros). With `init_with_zeros=false` the caller must have zeroed the array
    and every trace must have the same zero pattern, otherwise values written by one trace would be left in
    positions another trace treats as zero.
"""
function trace_function_expr(guards, J::AbstractMatrix{Node}, vars; in_place::Bool=false, init_with_zeros::Bool=true)
    body, nguard_statements, arguments = trace_function_body(guards, J, vars; in_place=in_place, init_with_zeros=init_with_zeros)
    expr = :(($(arguments...),) -> begin
        @inbounds begin
            $body
        end
    end)
    return expr, length(body.args), nguard_statements
end

# The statements of the generated function: guards, then the Jacobian. Returns the body, the number of statements that
# belong to the guards, and the names of the function's arguments.
function trace_function_body(guards, J::AbstractMatrix{Node}, vars; in_place::Bool=false, init_with_zeros::Bool=true)
    input_names, node_to_index = FD.variable_names((vars,))
    node_to_var = IdDict{Node,Union{Symbol,Real,Expr}}()
    body = Expr(:block)
    failure = in_place ? :(return false) : :(return nothing)
    for (cond, expected) in guards
        b, v = FD.function_body!(cond, node_to_index, node_to_var)
        append!(body.args, b.args)
        push!(body.args, expected ? :($v || $failure) : :($v && $failure))
    end
    if in_place
        if init_with_zeros && any(FD.is_identically_zero, J)
            push!(body.args, :(fill!(result, 0)))
        end
    else
        push!(body.args, :(result = zeros(promote_type(Float64, eltype($(input_names[1]))), $(size(J)))))
    end
    nguard_statements = length(body.args) # statements before the first Jacobian entry: the guards (and a fill of the result)
    for (i, node) in enumerate(J)
        FD.is_identically_zero(node) && continue
        b, v = FD.function_body!(node, node_to_index, node_to_var)
        append!(body.args, b.args)
        push!(body.args, :(result[$i] = $v))
    end
    push!(body.args, in_place ? :(return true) : :(return result))
    arguments = in_place ? (:result, input_names...) : Tuple(input_names)
    return body, nguard_statements, arguments
end

# ------------------------------------------------------------------ big traces: compile the statements in chunks

"""
A trace function split into chunks of straight-line code. Julia's compile time grows faster than linearly with the size
of one function, so a function of 100,000 statements compiles much more slowly than 50 functions of 2,000. Values that
are computed in one chunk and used in a later one are passed through the workspace vector `w`; everything else stays in
local variables. Only for `in_place=true` and `Float64` inputs.
"""
struct ChunkedFunction
    chunks::Vector{Any}
    w::Vector{Float64}
end
function (c::ChunkedFunction)(result, x)
    for f in c.chunks
        f(c.w, result, x) === true || return false
    end
    return true
end

const BOOL_OPS = (:<, :>, :<=, :>=, :(==), :isnan, :isinf, :isfinite, :signbit)

function symbols_in!(set::Set{Symbol}, e)
    if e isa Symbol
        push!(set, e)
    elseif e isa Expr
        for a in e.args
            symbols_in!(set, a)
        end
    end
    return set
end

function chunked_trace_function(guards, J::AbstractMatrix{Node}, vars; init_with_zeros::Bool=true, chunk_size::Int=2000)
    body, nguard, arguments = trace_function_body(guards, J, vars; in_place=true, init_with_zeros=init_with_zeros)
    stmts = body.args
    nst = length(stmts)
    chunk_of(i) = (i - 1) ÷ chunk_size + 1
    nchunks = chunk_of(nst)

    defchunk = Dict{Symbol,Int}()
    isbool = Set{Symbol}()
    for (i, s) in enumerate(stmts)
        if s isa Expr && s.head === :(=) && s.args[1] isa Symbol
            defchunk[s.args[1]] = chunk_of(i)
            rhs = s.args[2]
            rhs isa Expr && rhs.head === :call && rhs.args[1] in BOOL_OPS && push!(isbool, s.args[1])
        end
    end
    needs = [Set{Symbol}() for _ in 1:nchunks]
    exported = Set{Symbol}()
    for (i, s) in enumerate(stmts)
        c = chunk_of(i)
        for sym in symbols_in!(Set{Symbol}(), s)
            dc = get(defchunk, sym, 0)
            if dc != 0 && dc < c
                push!(needs[c], sym)
                push!(exported, sym)
            end
        end
    end
    slot = Dict{Symbol,Int}(sym => k for (k, sym) in enumerate(sort!(collect(exported))))

    functions = Any[]
    for c in 1:nchunks
        blk = Expr(:block)
        for sym in sort!(collect(needs[c]))
            push!(blk.args, isbool_load(sym, slot[sym], sym in isbool))
        end
        for i in ((c-1)*chunk_size+1):min(c * chunk_size, nst)
            s = stmts[i]
            if s === stmts[end] # return true: every chunk ends with its own
                continue
            end
            push!(blk.args, s)
            if s isa Expr && s.head === :(=) && s.args[1] isa Symbol && haskey(slot, s.args[1])
                sym = s.args[1]
                push!(blk.args, sym in isbool ? :(w[$(slot[sym])] = $sym ? 1.0 : 0.0) : :(w[$(slot[sym])] = $sym))
            end
        end
        push!(blk.args, :(return true))
        expr = :((w, $(arguments...),) -> begin
            @inbounds begin
                $blk
            end
        end)
        push!(functions, RGF.RuntimeGeneratedFunction(FD, FD, expr))
    end
    return ChunkedFunction(functions, zeros(length(slot))), nst, nguard
end
isbool_load(sym, k, boolean) = boolean ? :($sym = w[$k] != 0.0) : :($sym = w[$k])

"""Compile the generated expression with RuntimeGeneratedFunctions (the default way)."""
function make_trace_function(guards, J::AbstractMatrix{Node}, vars; in_place::Bool=false, init_with_zeros::Bool=true)
    expr, nstatements, nguard = trace_function_expr(guards, J, vars; in_place=in_place, init_with_zeros=init_with_zeros)
    return RGF.RuntimeGeneratedFunction(FD, FD, expr), nstatements, nguard
end

# Traces with more Jacobian nodes than this are compiled in chunks of CHUNK_SIZE statements (measured: a 154,619 statement
# trace compiles in 23 s in chunks of 10,000 against 170 s as one function, and its calls are about 1.4 times slower).
const CHUNK_THRESHOLD = Ref(30_000)
const CHUNK_SIZE = Ref(10_000)

# fn(x) for a chunked function: allocate the result, fill it, return it (or nothing if a guard failed)
struct AllocatingChunked{F}
    f::F
    sz::Tuple{Int,Int}
end
(a::AllocatingChunked)(x) = (r = zeros(a.sz); a.f(r, x) ? r : nothing)

function build_trace_function(tj, guards, J::AbstractMatrix{Node}, jacobian_nodes::Int)
    if jacobian_nodes > CHUNK_THRESHOLD[]
        cf, n, ng = chunked_trace_function(guards, J, tj.vars; init_with_zeros=tj.init_with_zeros, chunk_size=CHUNK_SIZE[])
        return (tj.in_place ? cf : AllocatingChunked(cf, size(J))), n, ng
    end
    return make_trace_function(guards, J, tj.vars; in_place=tj.in_place, init_with_zeros=tj.init_with_zeros)
end

# ------------------------------------------------------------------ traces and the cache

mutable struct TraceStats
    guards::Int
    trace_nodes::Int      # nodes in the recorded expression graph
    jacobian_nodes::Int   # nodes in the Jacobian expression graph
    statements::Int       # statements in the generated function
    run_s::Float64        # executing f on Tracked values (recording)
    diff_s::Float64       # FD jacobian
    codegen_s::Float64    # building the Expr
    first_call_s::Float64 # first call of the generated function, which includes Julia's compilation of it
    guard_statements::Int # statements of the generated function that belong to the guards, not to the Jacobian
    jacobian_method::Symbol # :factoring (FD.jacobian), :reverse (FD.reverse_AD, one pass per output) or :constant
end
TraceStats() = TraceStats(0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0, :none)

mutable struct Trace
    fn::Any # fn(x) -> Jacobian or nothing   (in_place=false)   |   fn(result, x) -> true or false   (in_place=true)
    stats::TraceStats
    bits::BitVector # the recorded decisions, used by guard_mode=:replay: exactly nbits comparison outcomes...
    ints::Vector{Int} # ...and nints integers
    nbits::Int
    nints::Int
    sig::UInt64     # signature of the decisions, to find the trace for a new input quickly
    J::Matrix{Node}                    # the Jacobian expressions, kept so that `fn` can be regenerated without guards...
    guards::Vector{Pair{Node,Bool}}    # ...or with them, when the function switches between compiled guards and replay
    with_guards::Bool                  # does `fn` contain the guards?
end

# Do the decisions in the live log equal the ones recorded for this trace? (Checked after the signatures match.)
function same_decisions(tr::Trace, l::DecisionLog)
    (tr.nbits == l.nbits && tr.nints == l.nints) || return false
    a = tr.bits.chunks
    b = l.bits.chunks
    nfull = l.nbits >> 6
    @inbounds for i in 1:nfull
        a[i] == b[i] || return false
    end
    r = l.nbits & 63
    if r != 0
        mask = (UInt64(1) << r) - 1
        @inbounds (a[nfull+1] & mask) == (b[nfull+1] & mask) || return false
    end
    @inbounds for i in 1:l.nints
        tr.ints[i] == l.ints[i] || return false
    end
    return true
end

"""
    TracingJacobian(f; wrt=[], max_traces=64, in_place=false, init_with_zeros=true)

Jacobian of `f` computed by tracing. `in_place` is fixed when the traces are generated:
  - `in_place=false`: `tj(args...)` returns a new Jacobian matrix.
  - `in_place=true`:  `tj(result, args...)` fills the matrix `result`; `eval_jacobian!(tj, result, x)` does the same
    for a flat input vector `x` and allocates nothing when a cached trace is used.
"""
mutable struct TracingJacobian
    f::Any
    wrt::Vector{Int}     # which arguments to differentiate with respect to (all if empty)
    max_traces::Int
    vars::Vector{Node}   # one variable per scalar input
    wrt_cols::Vector{Int}
    traces::Vector{Trace} # most recently used first
    hits::Int
    misses::Int
    evictions::Int
    in_place::Bool
    init_with_zeros::Bool
    sizes::Vector{Any}            # size of each argument, nothing for a number
    ranges::Vector{UnitRange{Int}} # where each argument sits in the flat input vector
    jsize::Union{Nothing,Tuple{Int,Int}}
    guard_mode::Symbol            # :compiled, :replay or :auto (see below)
    replayer::Any                 # a Replayer, made at the first call (its type depends on the arguments)
    active::Symbol                # what is in use now: :compiled (the guard conditions are compiled into the trace) or :replay (f is run on Probe numbers)
    switch_info::String           # why :auto chose what it chose
end
"""
guard_mode:
  - `:compiled` the guard conditions are compiled into each trace, which are tried one after another.
  - `:replay`    f is run once on value-only numbers that log their decisions, and the trace with those decisions is looked up.
  - `:auto`      starts as `:compiled`. After each new trace it measures a call of the compiled trace, a replay lookup, and a
                 rejection of an older trace, and switches (once) to `:replay` if the estimated cost is lower.
"""
TracingJacobian(f; wrt=Int[], max_traces=64, in_place=false, init_with_zeros=true, guard_mode=:auto) =
    TracingJacobian(f, collect(Int, wrt), max_traces, Node[], Int[], Trace[], 0, 0, 0, in_place, init_with_zeros, Any[], UnitRange{Int}[], nothing,
        guard_mode in (:compiled, :replay, :auto) ? guard_mode : error("guard_mode must be :compiled, :replay or :auto"), nothing,
        guard_mode === :replay ? :replay : :compiled, "")

summary_line(tj::TracingJacobian) = "mode = $(tj.active), traces cached = $(length(tj.traces)), hits = $(tj.hits), misses (recordings) = $(tj.misses), evictions = $(tj.evictions)"

function trace_report(tr::Trace)
    s = tr.stats
    return "guards = $(s.guards), trace nodes = $(s.trace_nodes), Jacobian nodes = $(s.jacobian_nodes) [$(s.jacobian_method)], statements = $(s.statements); " *
           "seconds: run f = $(round(s.run_s, digits=3)), differentiate = $(round(s.diff_s, digits=3)), build code = $(round(s.codegen_s, digits=3)), compile + first call = $(round(s.first_call_s, digits=3))"
end

# FD.jacobian (graph factoring) on a recorded trace gave 3 to 10 times more operations than FD's own reverse accumulation
# (FD.reverse_AD) on the Jacobi SVD traces, and takes minutes on large ones. So: always compute the reverse mode rows
# (one pass per output, they share the forward nodes), and also the factored Jacobian when the trace is small enough
# for that to be affordable, and keep whichever has fewer operations.
const FACTORING_MAX_TRACE_NODES = Ref(20_000)
const MAX_STATEMENTS = Ref(typemax(Int)) # refuse to compile a Jacobian with more nodes than this (set by test runners)

operations(x) = count(FD.is_tree, FD.all_nodes(x))

function best_jacobian(outputs::Vector{Node}, wrtvars::Vector{Node}, trace_nodes::Int)
    nout, nwrt = length(outputs), length(wrtvars)
    all(FD.is_constant, outputs) && return fill(Node(0.0), nout, nwrt), :constant
    Jr = Matrix{Node}(undef, nout, nwrt)
    for (i, o) in enumerate(outputs)
        Jr[i, :] = FD.is_constant(o) ? zeros(Node, nwrt) : FD.reverse_AD(o, wrtvars)
    end
    trace_nodes > FACTORING_MAX_TRACE_NODES[] && return Jr, :reverse
    Jf = FD.jacobian(outputs, wrtvars)
    return operations(Jf) <= operations(Jr) ? (Jf, :factoring) : (Jr, :reverse)
end

function init!(tj::TracingJacobian, args)
    ranges = arg_ranges(args)
    tj.ranges = ranges
    tj.sizes = Any[a isa Number ? nothing : size(a) for a in args]
    tj.replayer = make_replayer(tj.f, args)
    n = last(ranges[end])
    tj.vars = [Node(Symbol("x", i)) for i in 1:n]
    which = isempty(tj.wrt) ? collect(1:length(args)) : tj.wrt
    tj.wrt_cols = reduce(vcat, [collect(ranges[i]) for i in which])
end

function record!(tj::TracingJacobian, args, x)
    st = TraceStats()
    targs = tracked_args(args, arg_ranges(args), x, tj.vars)
    rec = Recording()
    CURRENT[] = rec
    t0 = time_ns()
    out = try
        tj.f(targs...)
    finally
        CURRENT[] = nothing
    end
    st.run_s = (time_ns() - t0) / 1e9
    outputs = as_nodes(out)
    # in :compiled mode "guards" is the number of conditions compiled into the trace, in :replay mode it is the number of
    # recorded decisions (every comparison and integer conversion that was executed)
    st.guards = tj.active === :replay ? rec.log.nbits + rec.log.nints : length(rec.guards)
    st.trace_nodes = length(FD.all_nodes(outputs))

    t0 = time_ns()
    wrtvars = tj.vars[tj.wrt_cols]
    J, st.jacobian_method = best_jacobian(outputs, wrtvars, st.trace_nodes)
    st.diff_s = (time_ns() - t0) / 1e9
    st.jacobian_nodes = length(FD.all_nodes(J))
    st.jacobian_nodes > MAX_STATEMENTS[] && error("TRACE_TOO_BIG: the Jacobian has $(st.jacobian_nodes) nodes, more than MAX_STATEMENTS = $(MAX_STATEMENTS[])")

    if tj.in_place
        tj.jsize === nothing || tj.jsize == size(J) || error("different paths gave Jacobians of different sizes ($(tj.jsize) and $(size(J))), which an in place Jacobian cannot hold")
        tj.jsize = size(J)
    end

    t0 = time_ns()
    compiled_guards = tj.active === :compiled ? rec.guards : Pair{Node,Bool}[] # in :replay mode the trace has no guards
    fn, nstatements, nguard = build_trace_function(tj, compiled_guards, J, st.jacobian_nodes)
    st.codegen_s = (time_ns() - t0) / 1e9
    st.statements = nstatements
    st.guard_statements = nguard
    l = rec.log
    return Trace(fn, st, l.bits[1:l.nbits], l.ints[1:l.nints], l.nbits, l.nints, signature(l), J, rec.guards, tj.active === :compiled)
end

# A trace whose function contains its guards is used in replay mode: generate the function without them. This is only
# needed for traces that were recorded before the function switched to replay; it is compiled when first called.
function make_guard_free!(tj::TracingJacobian, tr::Trace)
    fn, _, _ = build_trace_function(tj, Pair{Node,Bool}[], tr.J, tr.stats.jacobian_nodes)
    tr.fn = fn
    tr.with_guards = false
    return nothing
end

# :auto mode. Called after a new trace was recorded and added to the cache. Measures, on the input the trace was recorded
# from,
#   t_fused   a call of the trace with its guards compiled in
#   t_probe   one replay: run f on Probe numbers and look the trace up by its decisions
#   t_reject  what it costs the older traces to reject this input (they stop at their first failing guard)
# and estimates the cost of a call in each mode:
#   compiled  t_fused + (k - 1)/2 * t_reject      the traces are tried one after another, the right one is on average half way
#   replay    t_probe + t_jacobian                f is run once, then the Jacobian code without the guards
# The time of the Jacobian code alone is estimated as the fused call times the fraction of its statements that are not
# guards. (An earlier version assumed one half; that was wrong for functions with few guards and a big Jacobian.)
# It switches (once) to :replay if that is cheaper. The averages include garbage collection, which matters because the
# replay allocates (whatever f allocates).
function maybe_switch!(tj::TracingJacobian, x, result)
    (tj.guard_mode === :auto && tj.active === :compiled) || return nothing
    tr = tj.traces[1]
    call(t) = tj.in_place ? t.fn(result, x) : t.fn(x)
    function best(thunk; n=128, rounds=3)
        thunk(); thunk()
        return sum(begin
            t0 = time_ns()
            for _ in 1:n
                thunk()
            end
            (time_ns() - t0) / n
        end for _ in 1:rounds) / rounds
    end
    t_fused = best(() -> call(tr))
    t_probe = best(() -> find_trace(tj, x))
    others = tj.traces[2:min(end, 4)]
    t_reject = isempty(others) ? 0.0 : best(() -> (for o in others
        try
            call(o)
        catch
        end
    end)) / length(others)
    k = length(tj.traces)
    compiled_cost = t_fused + (k - 1) / 2 * t_reject
    jacobian_fraction = 1 - tr.stats.guard_statements / max(tr.stats.statements, 1)
    replay_cost = t_probe + jacobian_fraction * t_fused
    tj.switch_info = @sprintf("%d traces: fused %.0f ns (%.0f%% of its statements are the Jacobian), replay lookup %.0f ns, rejection %.0f ns => estimated compiled %.0f ns vs replay %.0f ns",
        k, t_fused, 100 * jacobian_fraction, t_probe, t_reject, compiled_cost, replay_cost)
    replay_cost < compiled_cost && (tj.active = :replay)
    return nothing
end

# Replay mode: run f once on Probe numbers, which write their decisions to REPLAY_LOG, and find the cached trace with
# exactly those decisions. Returns its position, or nothing if there is none (or if f raises an error on this input, in
# which case recording the input will raise the same error).
function find_trace(tj::TracingJacobian, x)
    reset!(REPLAY_LOG)
    try
        run_replayer(tj.replayer, x) # one dynamic call, then everything inside is type stable
    catch e
        e isa InterruptException && rethrow()
        return nothing
    end
    sig = signature(REPLAY_LOG)
    for (k, tr) in enumerate(tj.traces)
        tr.sig == sig && same_decisions(tr, REPLAY_LOG) && return k
    end
    return nothing
end

# A trace was just recorded from x. Replaying f on x must give the recorded decisions; if not, f is not deterministic
# or reads a traced value without going through the overloaded operations.
function check_replay_agrees(tj::TracingJacobian, x)
    tj.active === :replay && find_trace(tj, x) === nothing &&
        error("replaying f on the input a trace was recorded from gave different decisions than the recording. Is f deterministic?")
    return nothing
end

# the flat input vector back into the arguments, which is only needed when a new trace has to be recorded
unflatten(tj::TracingJacobian, x) = Any[sz === nothing ? x[first(r)] : reshape(x[r], sz) for (sz, r) in zip(tj.sizes, tj.ranges)]

function remember!(tj::TracingJacobian, tr::Trace)
    pushfirst!(tj.traces, tr)
    while length(tj.traces) > tj.max_traces
        pop!(tj.traces) # least recently used
        tj.evictions += 1
    end
end
touch!(tj::TracingJacobian, k, tr) = k > 1 && (deleteat!(tj.traces, k); pushfirst!(tj.traces, tr)) # keep the most recently used first

"""
    eval_jacobian!(tj, result, x)

Fill `result` with the Jacobian at the flat input vector `x` (`tj` must have been made with `in_place=true`). When a
cached trace's guards hold this allocates nothing. Otherwise a new trace is recorded, which does allocate.
"""
function eval_jacobian!(tj::TracingJacobian, result::AbstractMatrix, x::AbstractVector{Float64})
    tj.in_place || error("this TracingJacobian was made with in_place=false. Use tj(args...) or eval_jacobian(tj, x)")
    isempty(tj.vars) && error("call tj(result, args...) once first so that the shapes of the arguments are known")
    length(x) == length(tj.vars) || error("the arguments must have the same total length every call")
    tj.jsize === nothing || size(result) == tj.jsize || error("result has size $(size(result)) but the Jacobian is $(tj.jsize)")
    if tj.active === :replay
        k = find_trace(tj, x)
        if k !== nothing
            tr = tj.traces[k]
            tr.with_guards && make_guard_free!(tj, tr)
            tr.fn(result, x) # a replay trace has no guards, so this always fills result
            tj.hits += 1
            touch!(tj, k, tr)
            return result
        end
    else
        for (k, tr) in enumerate(tj.traces)
            ok = try
                tr.fn(result, x)
            catch e
                e isa InterruptException && rethrow()
                false # a guard expression that is only defined on another path (for example log of a negative number)
            end
            if ok === true
                tj.hits += 1
                touch!(tj, k, tr)
                return result
            end
        end
    end
    tj.misses += 1
    tr = record!(tj, unflatten(tj, x), x)
    size(result) == tj.jsize || error("result has size $(size(result)) but the Jacobian is $(tj.jsize)")
    t0 = time_ns()
    ok = tr.fn(result, x) # the first call compiles the generated function
    tr.stats.first_call_s = (time_ns() - t0) / 1e9
    ok === true || error("the guards of a trace failed on the input it was recorded from. Is the function deterministic?")
    remember!(tj, tr)
    check_replay_agrees(tj, x)
    maybe_switch!(tj, x, result)
    return result
end

"""
    eval_jacobian(tj, x)

Jacobian at the flat input vector `x` as a newly allocated matrix (`tj` must have been made with `in_place=false`).
"""
function eval_jacobian(tj::TracingJacobian, x::AbstractVector{Float64})
    tj.in_place && error("this TracingJacobian was made with in_place=true. Use eval_jacobian!(tj, result, x)")
    isempty(tj.vars) && error("call tj(args...) once first so that the shapes of the arguments are known")
    length(x) == length(tj.vars) || error("the arguments must have the same total length every call")
    if tj.active === :replay
        k = find_trace(tj, x)
        if k !== nothing
            tr = tj.traces[k]
            tr.with_guards && make_guard_free!(tj, tr)
            tj.hits += 1
            touch!(tj, k, tr)
            return tr.fn(x)
        end
    else
        for (k, tr) in enumerate(tj.traces)
            r = try
                tr.fn(x)
            catch e
                e isa InterruptException && rethrow()
                nothing # a guard expression that is only defined on another path
            end
            if r !== nothing
                tj.hits += 1
                touch!(tj, k, tr)
                return r
            end
        end
    end
    tj.misses += 1
    tr = record!(tj, unflatten(tj, x), x)
    t0 = time_ns()
    r = tr.fn(x) # the first call compiles the generated function
    tr.stats.first_call_s = (time_ns() - t0) / 1e9
    r === nothing && error("the guards of a trace failed on the input it was recorded from. Is the function deterministic?")
    remember!(tj, tr)
    check_replay_agrees(tj, x)
    maybe_switch!(tj, x, nothing)
    return r
end

# Convenience: take the arguments themselves. `tj(args...)` or, if in_place, `tj(result, args...)`. Flattening the
# arguments allocates a vector, so use eval_jacobian!/eval_jacobian with a flat vector inside a loop.
function (tj::TracingJacobian)(args...)
    if tj.in_place
        result = args[1]
        rest = Base.tail(args)
        isempty(tj.vars) && init!(tj, rest)
        return eval_jacobian!(tj, result, flatten(rest))
    else
        isempty(tj.vars) && init!(tj, args)
        return eval_jacobian(tj, flatten(args))
    end
end

# ------------------------------------------------------------------ loops: chain the Jacobians of the body trace

"""
    loop_jacobian(body_tj, body, cond, s0, p; max_iter)

Differentiate a loop `while cond(s, p); s = body(s, p); end` with a runtime trip count. Only the *body* is traced
(`body_tj = TracingJacobian(body)`), so the code size does not depend on the number of iterations. At each
iteration the Jacobian of the body trace `[ds'/ds | ds'/dp]` is multiplied into the accumulated derivative
`D = d s / d (s0, p)`:  D <- (ds'/ds) D + [0 | ds'/dp].
Returns the final state, `D` and the number of iterations.
"""
function loop_jacobian(body_tj::TracingJacobian, body, cond, s0, p; max_iter=10^6)
    m, q = length(s0), length(p)
    s = collect(Float64, s0)
    pp = collect(Float64, p)
    D = [Matrix{Float64}(I, m, m) zeros(m, q)]
    iters = 0
    while cond(s, pp)
        iters += 1
        iters > max_iter && error("more than $max_iter iterations")
        J = body_tj(s, pp)
        snext = body(s, pp)
        D = J[:, 1:m] * D
        D[:, (m+1):end] .+= J[:, (m+1):end]
        s = snext
    end
    return s, D, iters
end

end # module

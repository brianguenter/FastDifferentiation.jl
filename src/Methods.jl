#This copyright notice only applies to this file.

# Copyright (c) 2020: Shashi Gowda, Yingbo Ma, Mason Protter, Julia Computing.

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import SpecialFunctions: gamma, loggamma, erf, erfc, erfcinv, erfi, erfcx,
    dawson, digamma, trigamma, invdigamma, polygamma,
    airyai, airyaiprime, airybi, airybiprime, besselj0,
    besselj1, bessely0, bessely1, besselj, bessely, besseli,
    besselk, hankelh1, hankelh2, polygamma, beta, logbeta

const monadic = [deg2rad, rad2deg, asind, log1p, acsch,
    acos, asec, acosh, acsc, cscd, log, tand, log10, csch, asinh,
    abs2, cosh, sin, cos, atan, cospi, cbrt, acosd, acoth, acotd,
    asecd, exp, acot, sqrt, sind, sinpi, asech, log2, tan, exp10,
    sech, coth, asin, cotd, cosd, sinh, abs, csc, tanh, secd,
    atand, sec, acscd, cot, exp2, expm1, atanh, gamma,
    loggamma, erf, erfc, erfcinv, erfi, erfcx, dawson, digamma,
    trigamma, invdigamma, polygamma, airyai, airyaiprime, airybi,
    airybiprime, besselj0, besselj1, bessely0, bessely1, isfinite]

const diadic = [max, min, hypot, atan, mod, rem, copysign,
    besselj, bessely, besseli, besselk, hankelh1, hankelh2,
    polygamma, beta, logbeta]
const previously_declared_for = Set([])

const basic_monadic = [-, +]
const basic_diadic = [+, -, *, /, //, \, ^]

# TODO: keep domains tighter than this
function number_methods(T, rhs1, rhs2, options=nothing)
    exprs = []

    skip_basics = options !== nothing ? options == :skipbasics : false
    only_basics = options !== nothing ? options == :onlybasics : false
    skips = Meta.isexpr(options, [:vcat, :hcat, :vect]) ? Set(options.args) : []

    for f in (skip_basics ? diadic : only_basics ? basic_diadic : vcat(basic_diadic, diadic))
        nameof(f) in skips && continue
        for S in previously_declared_for
            push!(exprs, quote
                (f::$(typeof(f)))(a::$T, b::$S) = $rhs2
                (f::$(typeof(f)))(a::$S, b::$T) = $rhs2
            end)
        end

        # TODO: modularize and make another macro?
        expr = quote
            (f::$(typeof(f)))(a::$T, b::$T) = $rhs2
            (f::$(typeof(f)))(a::$T, b::Real) = $rhs2
            (f::$(typeof(f)))(a::Real, b::$T) = $rhs2
            # (f::$(typeof(f)))(a::$T, b::Number) = $rhs2
            # (f::$(typeof(f)))(a::Number, b::$T) = $rhs2
        end

        push!(exprs, expr)
    end

    for f in (skip_basics ? monadic : only_basics ? basic_monadic : vcat(basic_monadic, monadic))
        nameof(f) in skips && continue
        if f === isfinite
            push!(exprs, :((f::$(typeof(f)))(a::$T) = true))
        else
            push!(exprs, :((f::$(typeof(f)))(a::$T) = $rhs1))
        end
    end
    push!(exprs, :(push!($previously_declared_for, $T)))
    Expr(:block, exprs...)
end

macro number_methods(T, rhs1, rhs2, options=nothing)
    number_methods(T, rhs1, rhs2, options) |> esc
end





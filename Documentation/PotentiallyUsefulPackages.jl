module Graph
using GraphViz

# other visualization packages to consider
# https://github.com/JuliaGraphs/GraphPlot.jl
# https://github.com/JuliaGraphs/NetworkLayout.jl

#This example code is from https://github.com/JuliaGraphs/GraphViz.jl

GraphViz.load("mygraph.dot")
show(dot"""
 digraph graphname {
     a -> b -> c;
     b -> d;
 }
""")
end #module

#links
# https://github.com/JuliaSymbolics/Symbolics.jl

#docs 
# https://juliadiff.org/DiffRules.jl/stable/
#
# """diffrule from https://github.com/JuliaDiff/DiffRules.jl, TermInterface api from https://github.com/JuliaSymbolics/TermInterface.jl"""
# function diffrules(expr)
#     diff = DiffRules.diffrule(:Base,:cos,expr)
#     println("expression $expr derivative $diff")
#     println("operation $(operation(diff)) arguments $(arguments(diff))")
# end
# export diffrules

# """examples from https://github.com/JuliaSymbolics/SymbolicUtils.jl"""
# function symbolicutils()
#     SymbolicUtils.show_simplified[] = true

#     @syms x::Real y::Real z::Complex f(::Number)::Real

#     println(2x^2 - y + x^2) #simplification of 2x^2 + x^2 to 3x2
#     println(f(sin(x)^2 + cos(x)^2) + z) #simplification of sin(x)^2 + cos(x)^2 to 1
#     r = @rule sinh(im * ~x) => sin(~x)
#     println(r(sinh(im * y))) #should print sin(y)

#     println(simplify(cos(y)^2 + sinh(im*y)^2, RuleSet([r]))) #should print 1
# end
# export symbolicutils
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
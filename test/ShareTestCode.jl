module FDTests
using StaticArrays
using Memoize
using DataStructures

using FastDifferentiation



const FD = FastDifferentiation
export FD

include("TestPrograms/TestCode.jl")

#Each line in the postorder listing has two spaces at the end which causes a line break. Don't delete the trailing spaces or the formatting will get messed up.

end #module

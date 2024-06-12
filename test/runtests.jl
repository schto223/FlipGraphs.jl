using Test
using FlipGraphs

if isempty(ARGS)
    tests = [
        "polygonTriangulations.jl",
        "flipGraph_planar.jl",
        "deltaComplex.jl",
        "flipGraph.jl",
        "show.jl",
        "exporting.jl",
        "generalUtilities.jl"
    ]
else
    tests = ARGS
end

for test in tests
    include(test)
end
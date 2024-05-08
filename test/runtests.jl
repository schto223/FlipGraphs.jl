using Test
using FlipGraphs

if isempty(ARGS)
    tests = [
        "polygonTriangulations.jl",
        "flipGraph_planar.jl",
        "deltaComplex.jl"
    ]
else
    tests = ARGS
end

for test in tests
    include(test)
end
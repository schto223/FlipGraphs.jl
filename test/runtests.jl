using Test, Documenter
using FlipGraphs

if isempty(ARGS)
    tests = [
        "polygonTriangulations.jl",
        "flipGraphPlanar.jl",
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

doctest(FlipGraphs)

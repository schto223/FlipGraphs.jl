using Test
using FlipGraphs

if isempty(ARGS)
    tests = [
        "planarTriangulations.jl"
    ]
else
    tests = ARGS
end

for test in tests
    include(test)
end

using FlipGraphs, Random

function square(A)
    return A*A
end

g = triangulated_polygon(20)

A = adjacency_matrix(g)
B = Matrix{Int32}(A)

d = diameter(B)

#@time FlipGraphs.diameter(A)
#@time FlipGraphs.diameter(A)
#@time FlipGraphs.diameter(B)
#@time FlipGraphs.diameter(B)


using FlipGraphs




n = 3
g = triangulatedPolygon(n)

#drawPNG(g)

G = construct_FlipGraph(g, false)

drawPNG(G)

for n in 3:10
    g = triangulatedPolygon(n)
    G = construct_FlipGraph(g, false)
    drawPNG(G, string("flipG-",n))
    G = construct_FlipGraph(g, true)
    drawPNG(G, string("flipReduced-",n))
end
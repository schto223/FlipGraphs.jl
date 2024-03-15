
using FlipGraphs




n = 6
g = triangulatedPolygon(n)

#drawPNG(g)

G = construct_FlipGraph(g)

drawPNG(G)
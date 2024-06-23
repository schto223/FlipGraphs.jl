
using FlipGraphs, Random
using GraphPlot: gplot, cycle_graph
using Cairo#, KernelDensityEstimatePlotting

import Fontconfig
import Compose: px, PNG, draw

function drawPNG(g::TriangulatedPolygon, fileName::String ="triangulatedPolygon", drawLabels::Bool = false )
    n = g.n
    a = 2*π/n 
    x = [cos(π/2 -a + a*i) for i = 1:n]
    y = [-sin(π/2 -a + a*i) for i = 1:n]
    nodeLabel = 1:n
    if drawLabels
        draw(PNG("img/"*fileName * ".png", 500px, 500px), gplot(g, x, y, nodelabel=nodeLabel))
    else
        draw(PNG("img/"*fileName*".png", 500px, 500px), gplot(g,x,y))
    end
end

function drawPNG(G::FlipGraphPlanar, fileName::String ="flipGraph" , drawLabels::Bool=false)
    n = nv(G)#len(G.V)
    nodeLabel = 1:n
    if drawLabels
        draw(PNG("img/"*fileName*".png", 1000px, 1000px), gplot(G, nodelabel=nodeLabel))
    else
        draw(PNG("img/"*fileName*".png", 1000px, 1000px), gplot(G))
    end
end

"""
    drawPNG(G::FlipGraph, fileName::String = "flipGraph" , drawLabels::Bool = false)

Create a PNG image of the FlipGraph `G`.

"""
function drawPNG(G::FlipGraph, fileName::String ="flipGraph" , drawLabels::Bool=false)
    n = nv(G)
    nodeLabel = 1:n
    if drawLabels
        draw(PNG("img/"*fileName*".png", 1000px, 1000px), gplot(G, nodelabel=nodeLabel, layout=spectral_layout))
    else
        draw(PNG("img/"*fileName*".png", 1000px, 1000px), gplot(G))
    end
end

#G = flipgraph_modular(1,2)
#drawPNG(G, "G_0-3")

#for i in 3:15
#    export_gml("gml_exports/G-$i", flipgraph_planar(i), :diameter);
#end

for (g,p) in [(0,3),(0,4),(0,5), (0,6), (0,7), (0,8), (1,2), (1,4), (1,5), (2,2), (2,3), (1,6)]
    export_gml("gml_exports/closed_surfaces/G~$g-$p", flipgraph_modular(g,p,labeled_points=false), :diameter);
end
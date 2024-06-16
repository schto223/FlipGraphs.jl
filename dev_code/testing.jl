
using FlipGraphs, Random
using GraphPlot#: gplot, cycle_graph
#using Cairo, KernelDensityEstimatePlotting

#import Fontconfig
#import Compose: px, PNG, draw

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

g = triangulated_polygon(6)
G = flipgraph(g, modular=true)

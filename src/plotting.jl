using GraphPlot#: gplot, cycle_graph
import Cairo, Fontconfig
import Compose: px, PNG, draw

export drawPNG, plot

"""
    drawPNG(G::FlipGraphPlanar, fileName::String = "flipGraph" , drawLabels::Bool = false)

Create a PNG image of the FlipGraph `G`.

"""
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
    drawPNG(G::FlipGraphPlanar, fileName::String = "flipGraph" , drawLabels::Bool = false)

Create a PNG image of the FlipGraph `G`.

"""
function drawPNG(G::FlipGraph, fileName::String ="flipGraph" , drawLabels::Bool=false)
    n = nv(G)
    nodeLabel = 1:n
    if drawLabels
        draw(PNG("img/"*fileName*".png", 2000px, 2000px), gplot(G, nodelabel=nodeLabel, layout=spectral_layout))
    else
        draw(PNG("img/"*fileName*".png", 2000px, 2000px), gplot(G))
    end
end


"""
    drawPNG(g::TriangulatedPolygon, fileName::String = "triangulatedPolygon" , drawLabels::Bool = false)

Create a PNG image of the triangulated Polygon `g`.

"""
function drawPNG(g::TriangulatedPolygon, fileName::String ="triangulatedPolygon", drawLabels::Bool = false )
    n = g.n
    a = 2*π/n 
    x = [cos(π/2 -a + a*i) for i = 1:n]
    y = [-sin(π/2 -a + a*i) for i = 1:n]
    nodeLabel = 1:n
    if drawLabels
        draw(PNG("img/"*fileName * ".png", 1000px, 1000px), gplot(g, x, y, nodelabel=nodeLabel))
    else
        draw(PNG("img/"*fileName*".png", 1000px, 1000px), gplot(g,x,y))
    end
end

"""
    plot(g::TriangulatedPolygon, drawLabels::Bool = false)
"""
function plot(g::TriangulatedPolygon , drawLabels::Bool = false)
    n = nv(g)
    a = 2*π/n 
    x = [cos(π/2 -a + a*i) for i = 1:n]
    y = [-sin(π/2 -a + a*i) for i = 1:n]
    nodeLabel = collect(1:n)
    if drawLabels
        gplot(g, x, y, nodelabel=nodeLabel)
    else
        gplot(g, x, y)
    end
end


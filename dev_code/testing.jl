
using FlipGraphs, Random

#using Graphs
using Karnak
using NetworkLayout
using Colors

g = G
@drawsvg begin
    background("black")
    sethue("grey40")
    fontsize(8)
    drawgraph(g, 
        layout=stress, 
        vertexlabels = 1:nv(g),
        vertexfillcolors = 
            [RGB(rand(3)/2...) 
               for i in 1:nv(g)]
    )
end 600 400



#G = construct_FlipGraph_planar(6, true)
#for i in 1:nv(G)
#    drawPNG(G.V[i], string(i))
#    i += 1
#end


#HD = holeyDeltaComplex(1) 
#G = construct_FlipGraph(HD, 10, true)

#d = diameter_triangulation(D)

#for n in 3:15
#    g = triangulatedPolygon(n)
    #G = construct_FlipGraph(g, false)
    #drawPNG(G, string("flipG-",n))
#    G = construct_FlipGraph(g, true)
#    println(string(n,"  ", ne(G)))
    #drawPNG(G, string("flipReduced-",n))
#end


#D = createDeltaComplex(2,10)
#
#display(D)
#flip!(D,3)
#
#println()
#
#display(D)
#diameter(D)
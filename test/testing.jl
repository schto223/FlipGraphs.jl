
using FlipGraphs


for n in 3:15
    g = triangulatedPolygon(n)
    #G = construct_FlipGraph(g, false)
    #drawPNG(G, string("flipG-",n))
    G = construct_FlipGraph(g, true)
    println(string(n,"  ", ne(G)))
    #drawPNG(G, string("flipReduced-",n))
end


#D = createDeltaComplex(2,10)
#
#display(D)
#flip!(D,3)
#
#println()
#
#display(D)
#diameter(D)
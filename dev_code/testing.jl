using FlipGraphs, Random
using Colors
import Graphs: Edge

include("plotting.jl")

function walk_unique_nodes()
    G= flipgraph_planar(8)
    j = 68
    g = G.V[68]
    V = [68]    
    while true
        i = j
        #drawPNG(g,"g-8-$i"; nodecolor = colorant"black", inedgecolor = colorant"steelblue1", outedgecolor = colorant"grey10")
        j = findfirst(k -> !(k in V), neighbors(G,i))
        if !isnothing(j)
            push!(V, neighbors(G,i)[j])
            j = neighbors(G,i)[j]
        else
            break
        end
        println(j)
    end
    j=0
    for i in V
        j+=1
        drawPNG(G.V[i],"g-8-$j"; nodecolor = colorant"black", inedgecolor = colorant"steelblue1", outedgecolor = colorant"grey10")
    end
end

function dostuff()
    G = flipgraph_planar(8)
    g = G.V[68]
    Random.seed!(100)
    E = filter(e -> is_inner(g, e), edges(g))
    shuffle!(E)
    j=6
    while !isempty(E)
        drawPNG(g,"g-8-$j"; nodecolor = colorant"black", inedgecolor = colorant"steelblue1", outedgecolor = colorant"grey10")
        remove_edge!(g, popfirst!(E))
        j-=1
    end
    drawPNG(g,"g-8-$j"; nodecolor = colorant"black", inedgecolor = colorant"steelblue1", outedgecolor = colorant"grey10")
end

function dosuff2()
    g = flipgraph_planar(8).V[37]
    #g = flipgraph_planar(8).V[68]
    n = g.n
    m = ne(g)
    a = 2*π/n
    x = [cos(π/2 -a + a*i) for i = 1:n]
    y = [-sin(π/2 -a + a*i) for i = 1:n]
    nodefillc = [nodecolor = colorant"black" for i in 1:n]
    edgestrokec  = [( is_outer(g,e) ? colorant"grey10" : colorant"steelblue1" ) for e in edges(g)]
    #edgestrokec[7] = colorant"red"
    #edgestrokec[5] = colorant"orange"
    #edgestrokec[11] = colorant"orange"
    #edgestrokec[3] = colorant"orange"
    #edgestrokec[8] = colorant"orange"

    draw(PNG("img/gt.png", 1000px, 1000px), gplot(g,x,y, nodefillc=nodefillc, edgestrokec =edgestrokec, edgelinewidth=1.0, EDGELINEWIDTH=3.0 ))
end

function printallvertices()
    n = 8
    G= flipgraph_planar(n) 
    for i in eachindex(G.V)
        drawPNG(G.V[i],"g-$n-$i"; nodecolor = colorant"black", inedgecolor = colorant"steelblue1", outedgecolor = colorant"grey10")
    end
end


function export_withclasses(g,p)
    G = flipgraph_modular(g,p)
    Gtilde = flipgraph_modular(g,p, labeled_points=false)

    vats = [Dict() for i in 1:nv(G)]
    for i in eachindex(G.V)
        for j in eachindex(Gtilde.V)
            if is_isomorphic(get_vertex(G, i).D, get_vertex(Gtilde,j).D, labeled_points=false)
                vats[i]["equivalence_class"] = j
                break
            end
        end
    end
    export_gml(string("gml_exports/closed_surfaces/G_",g,"-",p), G; vertex_attributes=vats, add_diameter=true);  
end


function export_withclasses(p)
    G = flipgraph_planar(p)
    Gtilde = flipgraph_planar(p, modular=true)

    vats = [Dict() for i in 1:nv(G)]
    for i in eachindex(G.V)
        for j in eachindex(Gtilde.V)
            if is_isomorphic(get_vertex(G, i), get_vertex(Gtilde,j))
                vats[i]["equivalence_class"] = j
                break
            end
        end
    end
    export_gml(string("gml_exports/convex_polygon/G_",p), G; vertex_attributes=vats, add_diameter=true);  
end


#printallvertices()

#println(string("start at ", (time()-3600*24*365*54-3600*24*189-3600*20)/60))
#export_gml("gml_exports/convex_polygon/G~17", flipgraph_planar(17,modular=true), add_diameter=true);
#println("done")


for i in 13:15
    #export_withclasses(i)
    export_gml("gml_exports/convex_polygon/G-$i", flipgraph_planar(i), add_diameter=true);
end

for i in 3:15
    export_gml("gml_exports/convex_polygon/G~$i", flipgraph_planar(i,modular=true), add_diameter=true);
end

for (g,p) in [(0,3),(0,4),(0,5), (0,6), (0,7), (0,8), (1,2), (1,3), (1,4), (1,5), (2,2), (2,3), (1,6)]
    export_gml("gml_exports/closed_surfaces/G~$g-$p", flipgraph_modular(g,p,labeled_points=false), add_diameter=true);
end

for (g,p) in [(0,3),(0,4),(0,5), (0,6), (1,1), (1,2), (1,4), (2,1), (2,2), (2,3)]
    export_withclasses(g,p)
    #export_gml("gml_exports/closed_surfaces/G_$g-$p", flipgraph_modular(g,p,labeled_points=true), add_diameter=true);
end




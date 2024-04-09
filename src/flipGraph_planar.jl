

export FlipGraph, drawPNG, plot, construct_FlipGraph

"""
    struct FlipGraph <: AbstractGraph{Int}
A Graph representing the FlipGraph of a convex polygon.
Vertices are different triangulations of the same convex polygon.
Two vertices are linked by an edge, if the respective graphs differ only by a single flip.

"""
struct FlipGraph <: AbstractGraph{Int}    
    V::Array{TriGraph,1}
    adjList::Array{Array{Int,1},1}

    function FlipGraph(g::TriGraph)
        new(Array{TriGraph, 1}([g]), Array{Array{Int,1},1}([[]]))
    end

    function FlipGraph()
        new(Array{TriGraph,1}(), Array{Array{Int,1},1}())
    end

end


function edges(G::FlipGraph)
    E = collect(Edge(i,j) for i = 1:len(G.V) for j in G.adjList[i])
    return filter!(e->(src(e)>dst(e)), E)
end 

edgetype(G::FlipGraph) = SimpleEdge{Int}
has_edge(G::FlipGraph, e::Edge) = (dst(e) ∈ G.adjList[src(e)])
has_edge(G::FlipGraph, s, d) = (d ∈ G.adjList[s])
has_vertex(G::FlipGraph, v) = (1 <= v && v <= nv(G))
inneighbors(G::FlipGraph, v) = G.adjList[v]
ne(G::FlipGraph) = sum(size(G.adjList[i],1) for i=1:len(G.adjList))÷2
nv(G::FlipGraph) = len(G.V)
outneighbors(G::FlipGraph,v) = G.adjList[v]
vertices(G::FlipGraph) = G.V
is_directed(G::FlipGraph) = false
is_directed(::Type{FlipGraph}) = false


function add_edge!(G::FlipGraph, v, w) 
    if !has_edge(G, v, w) && v!=w
        push!(G.adjList[v],w)
        push!(G.adjList[w],v)
    end
end

function add_vertex!(G::FlipGraph, g::TriGraph) 
    push!(G.V,g)
    push!(G.adjList,[])
end

function remove_edge!(G::FlipGraph, e::Edge)
    deleteat!(G.adjList[src(e)], findfirst(x->x==dst(e), G.adjList[src(e)]))
    deleteat!(G.adjList[dst(e)], findfirst(x->x==src(e), G.adjList[dst(e)]))
end


"""
    drawPNG(G::FlipGraph, fName::String ="flipGraph" , drawLabels::Bool=false)

Creates a PNG image of the Graph G.

"""
function drawPNG(G::FlipGraph, fName::String ="flipGraph" , drawLabels::Bool=false)
    n = nv(G)#len(G.V)
    nodeLabel = 1:n
    if drawLabels
        draw(PNG("img/"*fName*".png", 1000px, 1000px), gplot(G, nodelabel=nodeLabel))
    else
        draw(PNG("img/"*fName*".png", 1000px, 1000px), gplot(G))
    end
end



"""
    function construct_FlipGraph(g::TriGraph, reduce::Bool=true)
    
        Returns the **FlipGraph** for the triangulated Polygon g.  

    If reduce is true, then vertices are the classes of isomorphisms up to renaming the vertices. Each class is represented by one of its elements.
    If reduce is false, then each vertex is a different triangulation of the initial graph g.

"""
function construct_FlipGraph(g::TriGraph, reduce::Bool=true)
    G = FlipGraph()
    if reduce
        sigpi = mcKay(g)[1]
        g = rename_vertices(g, sigpi)
    end
    add_vertex!(G, g)
    
    queue = Array{Tuple{TriGraph,Int},1}()
    push!(queue,(g,1))
    
    if reduce
        while !isempty(queue)
            g, ind_g = popfirst!(queue)
            for e in g.E
                if flippable(g,e)
                    gg = flip(g,e)
                    sigma_pis = mcKay(gg)
                    newGraph = true
                    for i = 1:len(G.V)
                        if is_isomorph(G.V[i], gg, sigma_pis)
                            add_edge!(G, ind_g, i)
                            newGraph = false
                            break
                        end
                    end
                    if newGraph
                        add_vertex!(G, rename_vertices(gg, sigma_pis[1]))
                        add_edge!(G, ind_g, len(G.V))
                        push!(queue, (G.V[end], len(G.V)))
                    end
                end
            end
        end
    else 
        while !isempty(queue)
            g, ind_g = popfirst!(queue)
            for e in g.E 
                if flippable(g,e)
                    gg = flip(g,e)
                    newGraph = true
                    for i = 1:len(G.V)
                        if all(issetequal(G.V[i].adjList[j], gg.adjList[j]) for j in 1:nv(g))
                            add_edge!(G, ind_g, i)
                            newGraph = false
                            break
                        end
                    end
                    if newGraph
                        add_vertex!(G, gg)
                        add_edge!(G, ind_g, len(G.V))
                        push!(queue, (G.V[end], len(G.V)))
                    end
                end
            end
        end
    end

    return G
end

function construct_FlipGraph(n::Int, reduce::Bool=true)
    return construct_FlipGraph(triangulatedPolygon(n), reduce)
end


function rename_vertices(g::TriGraph, sigma_pi::Array{Int,1})
    gg = TriGraph(g.n, Array{Edge{Int},1}(), Array{Array{Int,1},1}([[] for i in 1:g.n]))
    for e in g.E
        add_edge!(gg,sigma_pi[e.src], sigma_pi[e.dst])
    end
    return gg
end

"""

"""
function is_isomorph(g1::TriGraph, g2::TriGraph, sigma_pis::Array{Array{Int,1},1})
    if sort(degrees(g1)) != sort(degrees(g2))
        return false
    end
    for sig in sigma_pis
        if all(e-> has_edge(g1, sig[src(e)], sig[dst(e)]), g2.E)
            return true
        end
    end
    return false
end

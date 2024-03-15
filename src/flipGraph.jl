

export FlipGraph, drawPNG, plot, construct_FlipGraph

struct FlipGraph <: AbstractGraph{Integer}    
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
    E = collect(Edge(i,j) for i = 1:len(G.V) for j in G.adjList)
    return filter!(e->(src(e)>dst(e)), E)
end 

edgetype(G::FlipGraph) = SimpleEdge{Int}
has_edge(G::FlipGraph, e::Edge) = (dst(e) ∈ G.adjList[src(e)])
has_edge(G::FlipGraph, s, d) = (d ∈ G.adjList[s])
has_vertex(G::FlipGraph, v) = (1 <= v && v <= nv(G))
inneighbors(G::FlipGraph, v) = G.adjList[v]
ne(G::FlipGraph) = size(sum(size(G.adjList[i],1) for i=1:G.adjList) ,1)/2
nv(G::FlipGraph) = len(G.V)
outneighbors(G::FlipGraph,v) = G.adjList[v]
vertices(G::FlipGraph) = G.V
is_directed(G::FlipGraph) = false
is_directed(::Type{FlipGraph}) = false


function add_edge!(G::FlipGraph, v, w) 
    push!(G.adjList[v],w)
    push!(G.adjList[w],v)
end

function add_vertex!(G::FlipGraph, g::TriGraph) 
    push!(G.V,g)
    push!(G.adjList,[])
end

function remove_edge!(G::FlipGraph, e::Edge)
    deleteat!(G.adjList[src(e)], findfirst(x->x==dst(e), G.adjList[src(e)]))
    deleteat!(G.adjList[dst(e)], findfirst(x->x==src(e), G.adjList[dst(e)]))
end



function drawPNG(G::FlipGraph, fName::String ="flipGraph" )
    n = len(G.V)
    nodeLabel = 1:n
    draw(PNG("img/"*fName*".png", 1000px, 1000px), gplot(G, nodelabel=nodeLabel))
end

function plot(G::FlipGraph)
    n= len(G.V)
    a = 2*π/n 
    x = [cos(π/2 -a + a*i) for i = 1:n]
    y = [-sin(π/2 -a + a*i) for i = 1:n]
    nodeLabel = 1:n
    gplot(G, x,y, nodelabel=nodeLabel)
end

function construct_FlipGraph(g::TriGraph)
    G = FlipGraph()
    sigpi = mcKay(g)[1]
    g = rename_vertices(g, sigpi)

    add_vertex!(G,g)
    

    queue = Array{Tuple{TriGraph,Int},1}()
    push!(queue,(g,1))
    
    while !isempty(queue)
        g, ind_g = popfirst!(queue)
        #ind_g = getindex(G.V, g)
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

    return G
end


function rename_vertices(g::TriGraph, sigma_pi::Array{Int,1})
    gg = TriGraph(g.n, Array{Edge{Int},1}(), Array{Array{Int,1},1}([[] for i in 1:g.n]))
    for e in g.E
        add_edge!(gg,sigma_pi[e.src], sigma_pi[e.dst])
    end
    return gg
end


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




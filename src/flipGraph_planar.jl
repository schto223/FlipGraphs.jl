

export FlipGraphPlanar, construct_FlipGraph

"""
    struct FlipGraphPlanar <: AbstractGraph{Int}
A Graph representing the FlipGraph of a convex polygon.
Vertices are different triangulations of the same convex polygon.
Two vertices are linked by an edge, if the respective graphs differ only by a single flip.
"""
struct FlipGraphPlanar <: AbstractGraph{Int}    
    V::Vector{TriangulatedPolygon}
    adjList::Vector{Vector{Int}}

    function FlipGraphPlanar()
        new(Array{TriangulatedPolygon,1}(), Vector{Vector{Int}}())
    end
end

function Base.show(io::IO, mime::MIME"text/plain", G::FlipGraphPlanar)
    print(io, string("FlipGraphPlanar with ", nv(G) , " vertices and ", ne(G), " edges")); 
end

function edges(G::FlipGraphPlanar)
    E = collect(Edge(i,j) for i = 1:len(G.V) for j in G.adjList[i])
    return filter!(e->(src(e)>dst(e)), E)
end 

edgetype(G::FlipGraphPlanar) = SimpleEdge{Int}
has_edge(G::FlipGraphPlanar, e::Edge) = (dst(e) ∈ G.adjList[src(e)])
has_edge(G::FlipGraphPlanar, s, d) = (d ∈ G.adjList[s])
has_vertex(G::FlipGraphPlanar, v) = (1 <= v && v <= nv(G))
inneighbors(G::FlipGraphPlanar, v) = G.adjList[v]
ne(G::FlipGraphPlanar) = sum(size(G.adjList[i],1) for i=1:len(G.adjList))÷2
nv(G::FlipGraphPlanar) = len(G.V)
outneighbors(G::FlipGraphPlanar,v) = G.adjList[v]
vertices(G::FlipGraphPlanar) = G.V
is_directed(G::FlipGraphPlanar) = false
is_directed(::Type{FlipGraphPlanar}) = false


function add_edge!(G::FlipGraphPlanar, v, w) 
    if !has_edge(G, v, w) && v!=w
        push!(G.adjList[v],w)
        push!(G.adjList[w],v)
    end
end

function add_vertex!(G::FlipGraphPlanar, g::TriangulatedPolygon) 
    push!(G.V,g)
    push!(G.adjList,[])
end

function remove_edge!(G::FlipGraphPlanar, e::Edge)
    deleteat!(G.adjList[src(e)], findfirst(x->x==dst(e), G.adjList[src(e)]))
    deleteat!(G.adjList[dst(e)], findfirst(x->x==src(e), G.adjList[dst(e)]))
end




"""
    construct_FlipGraph(g::TriGraph [, modular::Bool])
    
Construct the **FlipGraph** for the triangulated Polygon `g`.  

If `modular` is true, then vertices of the FlipGraph are the classes of isomorphisms up to renaming the vertices. Each class is represented by one of its elements.\\
If `modular` is false, then each vertex is a different triangulation of the initial graph `g`.\\
By default, `modular` is set to `true`

"""
function construct_FlipGraph(g::TriangulatedPolygon, modular::Bool=true)
    G = FlipGraphPlanar()
    if modular
        sigpi = mcKay(g)[1]
        g = rename_vertices(g, sigpi)
    end
    add_vertex!(G, g)
    
    queue = Array{Tuple{TriangulatedPolygon,Int},1}()
    push!(queue,(g,1))
    
    if modular
        while !isempty(queue)
            g, ind_g = popfirst!(queue)
            for e in edges(g)
                if is_flippable(g,e)
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
            for e in edges(g) 
                if is_flippable(g,e)
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

"""
    construct_FlipGraph(n::Int [, modular::Bool])

Construct the FlipGraph of a convex n-gon. 

By default, the FlipGraph is reduced to its modular form. To get the complete FlipGraph, set `modular` to false

```julia
    construct_FlipGraph(5, false)
```
"""
function construct_FlipGraph(n::Int, modular::Bool=true)
    return construct_FlipGraph(triangulatedPolygon(n), modular)
end


"""
    rename_vertices(g::TriGraph, sigma_pi::Array{Int,1})

Rename the vertices of `g` by applying the permutation `sigma_pi`.
"""
function rename_vertices(g::TriangulatedPolygon, sigma_pi::Vector{<:Integer})
    gg = TriangulatedPolygon(g.n, Vector{Vector{Int}}([[] for i in 1:g.n]))
    for e in edges(g)
        add_edge!(gg,sigma_pi[e.src], sigma_pi[e.dst])
    end
    return gg
end

"""
    is_isomorph(g1::TriGraph, g2::TriGraph, sigma_pis::Array{Array{Int,1},1})

Return true if g1 is identical to g2 up to a renaming of the vertices of g1 by one of the permutations in sigma_pis.

sigma_pis is 
"""
function is_isomorph(g1::TriangulatedPolygon, g2::TriangulatedPolygon, sigma_pis::Array{Array{Int,1},1})
    if sort(degrees(g1)) != sort(degrees(g2))
        return false
    end
    for sig in sigma_pis
        if all(e-> has_edge(g1, sig[src(e)], sig[dst(e)]), edges(g2))
            return true
        end
    end
    return false
end


function diameter(G::FlipGraphPlanar)
    return diameter(adjacency_matrix(G.adjList))
end
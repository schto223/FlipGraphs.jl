#export mcKay, relative_degrees

struct FGPVertex
    g :: TriangulatedPolygon
    degrees ::Vector{Int} 
    num_point_perms :: Int
    ind :: Int

    function FGPVertex(g::TriangulatedPolygon, ind::Integer)
        return new(g, degrees(g), 0, ind)
    end

    function FGPVertex(g::TriangulatedPolygon, ind::Integer, perms::Vector{Vector{T}}) where T<:Integer
        g = rename_vertices(g, perms[1])
        new(g, degrees(g), length(perms), ind)
    end
end


"""
    struct FlipGraphPlanar <: AbstractGraph{Int32}

A Graph representing the FlipGraph of a convex polygon. 

Vertices are different triangulations of the same convex polygon.

Two vertices are linked by an edge, if the respective graphs differ only by a single flip.
"""
struct FlipGraphPlanar <: AbstractGraph{Int32}
    V ::Vector{FGPVertex}
    adjList ::Vector{Vector{Int32}}
    modular ::Bool

    function FlipGraphPlanar(modular::Bool=false)
        new(Vector{FGPVertex}(), Vector{Vector{Int32}}(), modular)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", G::FlipGraphPlanar)
    print(io, string("FlipGraphPlanar with ", nv(G) , " vertices and ", ne(G), " edges")); 
end

"""
    edges(G::FlipGraphPlanar) ::Vector{Edge}

Construct an array containing all the edges in `G`.
"""
function edges(G::FlipGraphPlanar) ::Vector{Edge}
    E = collect(Edge(Int32(i),j) for i in eachindex(G.V) for j in G.adjList[i])
    return filter!(e -> (src(e) > dst(e)), E)
end 

edgetype(G::FlipGraphPlanar) = SimpleEdge{Int32}

"""
    has_edge(G::FlipGraphPlanar, e::Edge)

Return `true` if `e` is an edge in `G`.
"""
has_edge(G::FlipGraphPlanar, e::Edge) = (dst(e) ∈ G.adjList[src(e)])

"""
    has_edge(G::FlipGraphPlanar, s, d)

Return `true` if there is an edge between `s` and `d` in `G`.
"""
has_edge(G::FlipGraphPlanar, s, d) = (d ∈ G.adjList[s])

"""
    has_vertex(G::FlipGraphPlanar, v)

Return `true` if `v` is a valid index of a vertex in `G`.
"""
has_vertex(G::FlipGraphPlanar, i::Integer) = (1 <= i <= nv(G))

"""
    has_vertex(G::FlipGraphPlanar, v::FGPVertex)

Return `true` if `v` is a vertex in `G`.
"""
has_vertex(G::FlipGraphPlanar, v::FGPVertex) = (v in G.V)


"""
    neighbors(G::FlipGraphPlanar, v::Integer) -> Vector{Int32}

Return a list of all the indices of vertices in `G`, that are adjacent to `v`.
"""
neighbors(G::FlipGraphPlanar, v::Integer) = G.adjList[v] :: Vector{Int32}
inneighbors(G::FlipGraphPlanar, v) = G.adjList[v]
outneighbors(G::FlipGraphPlanar,v) = G.adjList[v]

"""
    ne(G::FlipGraphPlanar)

Return the number of edges in `G`.
"""
ne(G::FlipGraphPlanar) = sum(size(G.adjList[i], 1) for i in eachindex(G.adjList))÷2

"""
    nv(G::FlipGraphPlanar)

Return the number of vertices in `G`.
"""
nv(G::FlipGraphPlanar) = length(G.V)

"""
    vertices(G::FlipGraphPlanar)

Return the List of all vertices in `G`.
"""
vertices(G::FlipGraphPlanar) = G.V


get_vertex(G::FlipGraphPlanar, i::Integer) = G.V[i]
is_directed(G::FlipGraphPlanar) = false
is_directed(::Type{FlipGraphPlanar}) = false

function add_edge!(G::FlipGraphPlanar, v, w) 
    if !has_edge(G, v, w) && v!=w
        push!(G.adjList[v],w)
        push!(G.adjList[w],v)
    end
end

function add_vertex!(G::FlipGraphPlanar, g::TriangulatedPolygon) 
    fgpv = FGPVertex(g, length(G.V)+1)
    push!(G.V, fgpv)
    push!(G.adjList,[])
    return fgpv
end

function add_vertex!(G::FlipGraphPlanar, g::TriangulatedPolygon, perms::Vector{Vector{T}}) where T<:Integer
    fgpv = FGPVertex(g, length(G.V)+1, perms)
    push!(G.V, fgpv)
    push!(G.adjList,[])
    return fgpv
end

#function remove_edge!(G::FlipGraphPlanar, e::Edge)
#    deleteat!(G.adjList[src(e)], findfirst(x->x==dst(e), G.adjList[src(e)]))
#    deleteat!(G.adjList[dst(e)], findfirst(x->x==src(e), G.adjList[dst(e)]))
#end


"""
    flipgraph(g::TriangulatedPolygon; kwargs..)
    
Construct the **FlipGraph** for the TriangulatedPolygon `g`.

# Arguments
- 'modular::Bool = false' : by default the whole flip graph is constructed. If modular is set to true, then only the modular flip graph is constructed.
In a modular flip graph, vertices of the FlipGraph are classes of isomorphisms up to renaming the vertices. 
Each class is then represented by one of its elements.
"""
function flipgraph(g::TriangulatedPolygon; modular::Bool = false)
    G = FlipGraphPlanar()
    if modular
        p = mcKay(g)
        fgpv = add_vertex!(G, g, p)
    else
        for i in 1:nv(g)
            sort!(g.adjList[i])
        end
        fgpv = add_vertex!(G, g)
    end
    
    queue = Vector{FGPVertex}()
    push!(queue, fgpv)
    
    if modular
        while !isempty(queue)
            fgpv = popfirst!(queue)
            g = fgpv.g
            for e in edges(g)
                if is_flippable(g,e)
                    gg = flip(g,e)
                    permutations = mcKay(gg)
                    num_perms = length(permutations)
                    degrees_sorted = sort(degrees(gg))
                    newGraph = true
                    for vg in G.V
                        if  vg.num_point_perms == num_perms && vg.degrees == degrees_sorted && is_isomorphic(vg, gg, permutations)
                            add_edge!(G, fgpv.ind, vg.ind)
                            newGraph = false
                            break
                        end
                    end
                    if newGraph
                        new_v = add_vertex!(G, gg, permutations)
                        add_edge!(G, fgpv.ind, new_v.ind)
                        push!(queue, new_v)
                    end
                end
            end
        end
    else 
        while !isempty(queue)
            fgpv = popfirst!(queue)
            g = fgpv.g
            for e in edges(g) 
                if is_flippable(g,e)
                    gg = flip(g,e)
                    newGraph = true
                    for j in 1:nv(gg)
                        sort!(gg.adjList[j])
                    end
                    degs = degrees(gg)
                    for v in G.V
                        if v.degrees == degs && all(v.g.adjList[j] == gg.adjList[j] for j in 1:nv(g))
                            add_edge!(G, fgpv.ind, v.ind)
                            newGraph = false
                            break
                        end
                    end
                    if newGraph
                        new_v = add_vertex!(G, gg)
                        add_edge!(G, fgpv.ind, new_v.ind)
                        push!(queue, new_v)
                    end
                end
            end
        end
    end
    return G
end

"""
    flipgraph_planar(n::Integer; modular=false)

Construct the `FlipGraph` of a convex `n`-gon. 

If `modular=true`, the FlipGraph is reduced to its modular form.

# Examples
```julia-repl
julia> flipgraph_planar(6)
FlipGraphPlanar with 14 vertices and 21 edges
```
"""
function flipgraph_planar(n::Integer; modular::Bool = false)
    return flipgraph(triangulated_polygon(n); modular = modular)
end

"""
    rename_vertices(g::TriangulatedPolygon, p::Vector{<:Integer})

Rename the vertices of `g` by applying the permutation `p`.
"""
function rename_vertices(g::TriangulatedPolygon, p::Vector{<:Integer})
    gg = TriangulatedPolygon(g.n)
    for i in 1:g.n
        gg.adjList[p[i]] = p[g.adjList[i]]
    end
    return gg
end


function is_isomorphic(g1::FGPVertex, g2::TriangulatedPolygon, permutations::Vector{Vector{T}}) where T<:Integer
    for p in permutations
        if all(i-> all(j-> has_edge(g1.g, p[i], p[j]), g2.adjList[i]) , eachindex(g2.adjList))
            return true
        end
    end
    return false
end

"""
    diameter(G::FlipGraphPlanar)

Compute the diameter of `G`.
"""
diameter(G::FlipGraphPlanar) = diameter(adjacency_matrix(G.adjList))

"""
    relative_degrees(g::TriangulatedPolygon, U::Vector{<:Integer}, V::Vector{<:Integer}) -> Vector{<:Integer}

Count for each vertex in `U`, the number of incident edges, which are also incident to an edge in `V`.
"""
function relative_degrees(g::TriangulatedPolygon, U::Vector{<:Integer}, V::Vector{<:Integer}) :: Vector{<:Integer}
    rdegs = zeros(Int32, length(U))
    for i in eachindex(U), j in V
        if has_edge(g, U[i], j)
            rdegs[i] += 1
        end
    end
    return rdegs
end


"""
    mcKay(g::TriangulatedPolygon) -> Vector{Vector{<:Integer}}

Apply *McKay's canonical graph labeling algorithm* in order to determine all possible permutations 
of the vertices which give a canonical isomorphism class representant.

Return a list of all possible canonical point relabeling permutations `p` such that the i-th point should be relabeled as the `p[i]`-th point
"""
function mcKay(g::TriangulatedPolygon) :: Vector{Vector{Int}}
    # renamed from sigma

    #split V into partitions according to their degrees from smallest to biggest
    function split(V::Vector{<:Integer}, degs::Vector{<:Integer}) 
        sV = Vector{Vector{Int}}()
        deg = 0 
        k = length(V)
        while k > 0
            W = Vector{Int}()
            for i in eachindex(degs)
                if degs[i] == deg
                    push!(W, V[i])
                    k -= 1
                end
            end
            if !isempty(W)
                push!(sV, W)
            end
            deg += 1
        end
        return sV
    end

    #replace partitions by partitions as long as there are 2 elements in the same partition ...
    #...that may be differentiated by their relative degrees to another partition
    function makeEquitable!(p::Vector{Vector{T}}, g::TriangulatedPolygon) where T<:Integer
        i = 1; j = 1
        while i <= length(p)
            rDegs = relative_degrees(g, p[i], p[j])
            if !all(x -> x==rDegs[1], rDegs) 
                newVs = split(p[i], rDegs)
                #replace the old partition by the new ones
                popat!(p,i)
                j = i
                for V in newVs
                    insert!(p, i, V)
                    i += 1
                end
                i = 1
            else 
                j += 1
                if j > length(p)
                    j = 1
                    i += 1
                end
            end
        end    
    end

    p = split(collect(1:g.n), degrees(g))
    makeEquitable!(p, g)
    if length(p) == g.n #there is only one canonical permutation
        return Vector{Vector{Int}}([invert_permutation(reduce(vcat, p))])
    end
    
    #split the first partition that has more than 2 elements 
    queue = Vector{Vector{Vector{Int}}}([p])
    leafs = Vector{Vector{Vector{Int}}}()
    while !isempty(queue)
        p = popfirst!(queue)
        i = 1
        while length(p[i]) == 1
            i += 1
        end
        for j in eachindex(p[i])
            pp = deepcopy(p)
            V = popat!(pp, i)
            insert!(pp, i, [popat!(V, j)])
            insert!(pp, i+1, V)
            makeEquitable!(pp, g)
            if length(pp) != g.n
                push!(queue, pp)
            else
                push!(leafs, pp)
            end
        end
    end

    return [invert_permutation(reduce(vcat, sigpi)) for sigpi in leafs]   #permutations
end
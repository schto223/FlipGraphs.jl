"""
    struct FlipGraphPlanar <: AbstractGraph{Int32}

A Graph representing the **flip graph** of a **convex polygon**. 

Vertices are different triangulations of the same convex polygon.
Two vertices are linked by an edge, if the respective graphs differ only by a single flip.
"""
struct FlipGraphPlanar <: AbstractGraph{Int32}
    V ::Vector{TriangulatedPolygon}
    adjList ::Vector{Vector{Int32}}
    modular ::Bool

    function FlipGraphPlanar(modular::Bool=false)
        new(Vector{TriangulatedPolygon}(), Vector{Vector{Int32}}(), modular)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", G::FlipGraphPlanar)
    print(io, string("FlipGraphPlanar with ", nv(G) , " vertices and ", ne(G), " edges")); 
end


"""
    edges(G::FlipGraphPlanar) ::Vector{Edge}

Construct an array containing all the edges in `G`.
"""
function edges(G::FlipGraphPlanar) :: Vector{Edge}
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
    has_vertex(G::FlipGraphPlanar, i::Integer) :: Bool

Return `true` if `i` is a valid index of a vertex in `G`.
"""
has_vertex(G::FlipGraphPlanar, i::Integer) ::Bool = (1 <= i <= nv(G))

"""
    has_vertex(G::FlipGraphPlanar, g::TriangulatedPolygon) :: Bool

Return `true` if `g` is a vertex in `G`. 

If `G` is a modular flip graph, this will only return `true` if `g` is the proper representant of the vertex.
"""
has_vertex(G::FlipGraphPlanar, g::TriangulatedPolygon) :: Bool = (g in G.V)


"""
    neighbors(G::FlipGraphPlanar, v::Integer) :: Vector{Int32}

Return a list of all the indices of vertices in `G`, that are adjacent to `v`.
"""
neighbors(G::FlipGraphPlanar, v::Integer) :: Vector{Int32} = G.adjList[v] 
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

"""
    get_vertex(G::FlipGraphPlanar, i::Integer) :: TriangulatedPolygon

Return the `i`-th vertex in the planar flip graph `G`.
"""
get_vertex(G::FlipGraphPlanar, i::Integer) :: TriangulatedPolygon = G.V[i]
is_directed(G::FlipGraphPlanar) = false
is_directed(::Type{FlipGraphPlanar}) = false

function add_edge!(G::FlipGraphPlanar, v, w) 
    if !has_edge(G, v, w) && v!=w
        push!(G.adjList[v],w)
        push!(G.adjList[w],v)
    end
end

function add_vertex!(G::FlipGraphPlanar, g::TriangulatedPolygon) 
    push!(G.V, g)
    push!(G.adjList,[])
end

"""
    degree(G::FlipGraphPlanar, v::Integer)

Return the **degree** of the `v`-th vertex in the graph `G`.
"""
function degree(G::FlipGraphPlanar, v::Integer) :: Int
    return length(G.adjList[v])
end

"""
    flipgraph(g::TriangulatedPolygon; kwargs..)
    
Construct the **`FlipGraph`** for the `TriangulatedPolygon` `g`.

# Arguments
- 'modular::Bool = false' : by default, the whole flip graph is constructed. If `modular` is set to `true`, then only the modular flip graph is constructed.
In a *modular flip graph*, vertices of the flip graph are classes of isomorphisms up to renaming the vertices. 
Each class is then represented by one of its elements.
"""
function flipgraph(g::TriangulatedPolygon; modular::Bool = false)
    if !modular
        return flipgraph_planar_labeledpoints(g)
    else
        return flipgraph_planar_unlabeledpoints(g)
    end
end

#We assume n>=4: The first two are always 2, hence they are skipped. we start at the third degree in the sorted list and only store the increases. Most will be 0.
#function degree_steps(g::TriangulatedPolygon{T}) :: Vector{T} where T<:Integer
#    d = sort!([length(g.adjList[i]) for i in 1:g.n])
#    return d[3:end] - d[2:end-1]
#end

"""
    flipgraph_planar_labeledpoints(g::TriangulatedPolygon)

Compute the **flip graph** with labeled points from the root vertex `g`.
"""
function flipgraph_planar_labeledpoints(g::TriangulatedPolygon)
    nvg = nv(g)
    G = FlipGraphPlanar()

    #D is a search tree which branches of for each degree of a vertex. The leafs are the ids of the vertices in the flip graph
    D = Vector{Any}(undef, nvg-2)
    add_vertex!(G, g)
    d = D
    degs = degrees(g)
    degs .-= 1
    for i in degs[1:end-2]
        d[i] = Vector{Any}(undef, nvg-2)
        d = d[i]
    end
    d[degs[end-1]] = 1

    queue::Vector{typeof(g)} = [g]
    numVG = 1
    v_id = 0
    while !isempty(queue)
        v_id += 1
        fgpv = popfirst!(queue)
        g = deepcopy(fgpv)
        for i in 1:nvg
            for j in fgpv.adjList[i]
                if i<j && i+1!=j && (i!=1 || j!=nvg)                        
                    i_new, j_new = flip_get_edge!(g,i,j)
                    newGraph = false
                    k = 1
                    d = D
                    deg = 0
                    while k < nvg
                        deg = length(g.adjList[k]) - 1
                        if isassigned(d, deg)
                            d = d[deg]
                        else
                            if k == nvg - 1 #it is a new vertex
                                numVG += 1
                                d[deg] = numVG
                                new_v = deepcopy(g)
                                add_vertex!(G, new_v)
                                add_edge!(G, v_id, numVG)
                                push!(queue, new_v)
                            else
                                d[deg] = Vector{Any}(undef, nvg-2)
                                d = d[deg]
                            end 
                            newGraph = true
                        end
                        k += 1
                    end
                    if !newGraph
                        add_edge!(G, v_id, d)
                    end
                    #revert the flip
                    flip!(g, i_new, j_new)
                end
            end
        end
    end
    return G
end

"""
    min_lexicographic_degrees(g::TriangulatedPolygon)

Return the lexicographically minimal `degrees(gi)` for the dihedral group of `g`.
"""
function min_lexicographic_degrees(g::TriangulatedPolygon)
    n = g.n
    degs = degrees(g)
    d2s = findall(x-> x==2, degs)
    best = vcat(d2s[1]:n, 1:d2s[1]-1)
    mirrored = true
    for i in eachindex(d2s)
        d = d2s[i]
        j = 1
        while j < n
            if mirrored  
                if degs[(d+n-1-j)%n + 1] != degs[best[j+1]]
                    if degs[(d+n-1-j)%n + 1] < degs[best[j+1]] #new best found
                        best = vcat(d2s[i]:-1:1, n:-1:d2s[i]+1)
                    end
                    mirrored = false
                    break
                end
            else
                if degs[(d+j-1)%n + 1] != degs[best[j+1]]
                    if degs[(d+j-1)%n + 1] < degs[best[j+1]]
                        best = vcat(d2s[i]:n, 1:d2s[i]-1)
                    end
                    j = 0
                    mirrored = true
                end
            end
            j += 1
            if j==n
                mirrored = false
            end
        end
    end
    return degs[best]
end

function flipgraph_planar_unlabeledpoints(g::TriangulatedPolygon)
    nvg = nv(g)
    G = FlipGraphPlanar()
    add_vertex!(G, g)

    #D is a search tree which branches of for each degree of a vertex. The leafs are the ids of the vertices in the flip graph
    d = D = Vector{Any}(undef, nvg-2)
    degs = min_lexicographic_degrees(g)
    degs .-= 1
    for i in degs[1:end-2]
        d[i] = Vector{Any}(undef, nvg-2)
        d = d[i]
    end
    d[degs[end-1]] = 1

    queue::Vector{typeof(g)} = [g]
    numVG = 1
    v_id = 0
    while !isempty(queue)
        v_id += 1
        fgpv = popfirst!(queue)
        g = deepcopy(fgpv)
        for i in 1:nvg
            for j in fgpv.adjList[i]
                if i<j && i+1!=j && (i!=1 || j!=nvg)                        
                    i_new, j_new = flip_get_edge!(g,i,j)
                    newGraph = false
                    k = 1
                    d = D
                    degs = min_lexicographic_degrees(g).-1
                    while k < nvg
                        deg = degs[k]
                        if isassigned(d, deg)
                            d = d[deg]
                        else
                            if k == nvg - 1 #it is a new vertex
                                numVG += 1
                                d[deg] = numVG
                                new_v = deepcopy(g)
                                add_vertex!(G, new_v)
                                add_edge!(G, v_id, numVG)
                                push!(queue, new_v)
                            else
                                d[deg] = Vector{Any}(undef, nvg-2)
                                d = d[deg]
                            end 
                            newGraph = true
                        end
                        k += 1
                    end
                    if !newGraph
                        add_edge!(G, v_id, d)
                    end
                    #revert the flip
                    flip!(g, i_new, j_new)
                end
            end
        end
    end
    return G
end

"""
    flipgraph_planar(n::Integer; modular=false) :: FlipGraphPlanar

Construct the `FlipGraphPlanar` of a convex `n`-gon. 

If `modular=true`, the flip graph is reduced to its modular form.

# Examples
```julia-repl
julia> flipgraph_planar(6)
FlipGraphPlanar with 14 vertices and 21 edges
```
"""
function flipgraph_planar(n::Integer; modular::Bool = false) :: FlipGraphPlanar
    return flipgraph(triangulated_polygon(Int8, n); modular = modular)
end

"""
    rename_vertices(g::TriangulatedPolygon, p::Vector{<:Integer}) :: TriangulatedPolygon

Return the `TriangulatedPolygon` obtained from renaming the vertices of the triangulated convex polygon `g` by applying the permutation `p`.
"""
function rename_vertices(g::TriangulatedPolygon, p::Vector{<:Integer}) :: TriangulatedPolygon
    gg = TriangulatedPolygon{eltype(g)}(g.n)
    for i in 1:g.n
        gg.adjList[p[i]] = p[g.adjList[i]]
    end
    return gg
end

"""
    rename_vertices!(g::TriangulatedPolygon, p::Vector{<:Integer}) :: TriangulatedPolygon

Rename the vertices of the triangulated convex polygon `g` by applying the permutation `p`.
"""
function rename_vertices!(g::TriangulatedPolygon, p::Vector{<:Integer})
    g.adjList[p] = g.adjList[:]
    for i in 1:g.n
        g.adjList[i] = p[g.adjList[i]]
    end
    return g
end

"""
    is_isomorphic(g1::TriangulatedPolygon, g2::TriangulatedPolygon, permutations::Vector{Vector{T}}) where T<:Integer

Check if `g2` is isomorphic to `g1` up to a relabeling of the vertices by one of the `permutations`.
"""
function is_isomorphic(g1::TriangulatedPolygon, g2::TriangulatedPolygon, permutations::Vector{Vector{T}}) where T<:Integer
    for p in permutations
        if all(i-> all(j -> has_edge(g1, p[i], p[j]), g2.adjList[i]), eachindex(g2.adjList))
            return true
        end
    end
    return false
end

"""
    is_isomorphic(g1::TriangulatedPolygon, g2::TriangulatedPolygon)

Check if `g1` and `g2` are isomorphic up to a relabeling of the vertices.
"""
function is_isomorphic(g1::TriangulatedPolygon, g2::TriangulatedPolygon)
    g = rename_vertices(g1, mcKay(g1, only_one=true)[1])
    permutations = mcKay(g2)
    for p in permutations
        if all(i-> all(j-> has_edge(g, p[i], p[j]), g2.adjList[i]) , eachindex(g2.adjList))
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
    mcKay(g::TriangulatedPolygon) :: Vector{Vector{<:Integer}}

Apply *McKay's canonical graph labeling algorithm* to determine all possible permutations 
of the vertices, which give a canonical isomorphism class representant.

Return a list of all possible canonical point-relabeling permutations `p` such that the i-th point should be relabeled as the `p[i]`-th point
"""
function mcKay(g::TriangulatedPolygon{T}; only_one::Bool =false) :: Vector{Vector{T}} where T<:Integer
    #split V into partitions according to their degrees from smallest to biggest
    function split(V::Vector{T}, degs::Vector{<:Integer}) :: Vector{Vector{T}} where T<:Integer
        sV = Vector{Vector{T}}()
        deg = 0 
        k = length(V)
        while k > 0
            W = Vector{T}()
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
            bo = true
            rd1 = relative_degree(g, p[i][1], p[j])
            for k in 2:length(p[i])
                if relative_degree(g, p[i][k], p[j]) != rd1
                    bo = false
                    break
                end
            end
            if !bo
                rDegs = relative_degrees(g, p[i], p[j])
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

    p = split(collect(T, 1:g.n), degrees(g))
    makeEquitable!(p, g)
    if length(p) == g.n #there is only one canonical permutation
        return Vector{Vector{T}}([invert_permutation(reduce(vcat, p))])
    end

    while only_one   
        i = 1
        while length(p[i]) == 1 #i = index of first partition that is not trivial
            i += 1 
        end 
        V = popat!(p, i)
        insert!(p, i, [popat!(V, 1)])
        insert!(p, i+1, V)
        makeEquitable!(p, g)
        if length(p) == g.n
            return Vector{Vector{T}}([invert_permutation(reduce(vcat, p))])
        end
    end
    
    #split the first partition that has more than 2 elements 
    queue = Vector{Vector{Vector{T}}}([p])
    leafs = Vector{Vector{T}}()
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
                push!(leafs, invert_permutation(reduce(vcat, pp)))
            end
        end
    end
    return leafs
end
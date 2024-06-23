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
    end
    nvg = nv(g)

    G = FlipGraphPlanar()
    p = mcKay(g, only_one=true)[1]
    add_vertex!(G, rename_vertices!(g,p))
    
    #D is a search tree which sorts vertices by their sorted degrees. 
    #Since arrays have to start at 1, and we do not want leading empty values, we look at the increase in the degree from vertex to vertex
    #max_steps gives an upper limit to the possible step size at that point.
    max_cum_degs = 2*nvg - 6 
    deg_steps = degree_steps(g)
    max_steps = zeros(Int, nvg-2)
    for i in eachindex(max_steps)
        max_steps[i] =  max_cum_degs ÷ (nvg-1-i)
    end
    D = Vector{Vector{Any}}(undef, 1 + max_cum_degs÷(nvg-2))
    d = D  
    for i in eachindex(deg_steps)
        if i == nvg-2
            d[deg_steps[i] + 1] = Int[1]
        else
            d[deg_steps[i] + 1] = Vector{Any}(undef, 1 + max_steps[i+1])
            d = d[deg_steps[i] + 1]
        end    
    end

    id  = 0
    depth  = 0
    numV  = 1
    depth_next_jump = 2
    depth_1  = 0
    depth_2 = 0
    depth_cutof = 0
    while id < length(G.V)
        id += 1
        if id == depth_next_jump
            depth += 1
            depth_cutof, depth_1, depth_2, depth_next_jump = depth_1, depth_2, depth_next_jump, numV
        end
        g_v = G.V[id]
        g = deepcopy(g_v)
        for i in 1:nvg
            for j in g_v.adjList[i]
                if  i<j && count(k in g.adjList[j] for k in g.adjList[i]) == 2 # i-j is a flippable edge  is_flippable(g ,i ,j)
                    i_new, j_new = flip_get_edge!(g,i,j)
                    permutations = mcKay(g, only_one=true)
                    d = D
                    newGraph = false
                    degs = sort!(degrees(g))
                    for k in eachindex(deg_steps)
                        deg_step = degs[k+2] - degs[k+1] + 1
                        if isassigned(d, deg_step)
                            d = d[deg_step]
                        else
                            if k == nvg - 2
                                d[deg_step] = Int[]
                                d = d[deg_step]
                            else
                                d[deg_step] = Vector{Vector{Any}}(undef, 1 + max_steps[k+1])
                                d = d[deg_step]
                            end
                            newGraph = true
                        end            
                    end
                    if !newGraph
                        newGraph = true
                        for v_id in d
                            if v_id >= depth_cutof
                                if is_isomorphic(G.V[v_id], g, permutations)
                                    add_edge!(G, v_id, id)
                                    newGraph = false
                                    break
                                end
                            end
                        end
                    end
                    if newGraph
                        new_v = rename_vertices(g, permutations[1])
                        add_vertex!(G, new_v)
                        numV += 1
                        add_edge!(G, id, numV)
                        push!(d, numV)
                    end
                    flip!(g, i_new, j_new)
                end
            end
        end
    end
    return G
end

#We assume n>=4: The first two are always 2, hence they are skipped. we start at the third degree in the sorted list and only store the increases. Most will be 0.
function degree_steps(g::TriangulatedPolygon{T}) :: Vector{T} where T<:Integer
    d = sort!([length(g.adjList[i]) for i in 1:g.n])
    return d[3:end] - d[2:end-1]
end

"""
    flipgraph_planar_labeledpoints(g::TriangulatedPolygon)

Compute the **flip graph** with labeled points from the root vertex `g`.
"""
function flipgraph_planar_labeledpoints(g::TriangulatedPolygon)
    nvg = nv(g)
    G = FlipGraphPlanar()
    push!(G.adjList, [])

    #D is a search tree which branches of for each degree of a vertex. The leafs are the ids of the vertices in the flip graph
    D = Vector{Any}(undef, nvg-2)
    push!(G.V, g)
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
                                push!(G.V, new_v)
                                push!(G.adjList, [])
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
    rename_vertices!(g1, mcKay(g1, only_one=true)[1])
    permutations = mcKay(g2)
    for p in permutations
        if all(i-> all(j-> has_edge(g1, p[i], p[j]), g2.adjList[i]) , eachindex(g2.adjList))
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
    relative_degrees(g::TriangulatedPolygon, U::Vector{<:Integer}, V::Vector{<:Integer}) :: Vector{<:Integer}

Count, for each vertex in `U`, the number of incident edges, which are also incident to an edge in `V`.
"""
function relative_degrees(g::TriangulatedPolygon, U::Vector{<:Integer}, V::Vector{<:Integer}) :: Vector{Int32}
    rdegs = zeros(Int32, length(U))
    for i in eachindex(U), j in V
        if has_edge(g, U[i], j)
            rdegs[i] += 1
        end
    end
    return rdegs
end

"""
    relative_degree(g::TriangulatedPolygon, u::Integer, V::Vector{<:Integer}) :: Vector{<:Integer}

Count the number of edges in `g` going from `u` to a vertex in `V`.
"""
function relative_degree(g::TriangulatedPolygon, u::Integer, V::Vector{T}) :: T where T<:Integer
    adj = g.adjList[u]
    return count(j->j in adj, V)
end

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
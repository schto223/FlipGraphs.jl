
export TriangulatedPolygon, triangulated_polygon
export edges, vertices, has_edge, has_vertex, neighbors, ne, nv, degrees
export flip, flip!, is_flippable

"""
    struct TriangulatedPolygon <: AbstractGraph{Integer} 

A structure representing a triangulation of a convex polygon.
"""
struct TriangulatedPolygon <: AbstractGraph{Integer}    
    n ::Int
    adjList ::Vector{Vector{Int}}

    function TriangulatedPolygon(n::Integer)
        new(n, Vector{Vector{Int}}([[] for i in 1:n]))
    end

    function TriangulatedPolygon(n::Integer, adjList::Vector{Vector{T}}) where {T<:Integer}
        new(n, adjList)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", g::TriangulatedPolygon)
    print(io, string("TriangulatedPolygon with ",g.n, " vertices, and adjacency list:")); println(io)
    print(io, g.adjList)
end

"""
    triangulated_polygon(n::Int) -> TriangulatedPolygon

Create a triangulated convex n-gon. 

Vertices are named from 1 to n in an anticlockwise manner.
The inside is triangulated in a zig-zag manner.
"""
function triangulated_polygon(n::Integer) :: TriangulatedPolygon
    n>=0 || throw(ArgumentError(string("The number of vertices(n) cannot be negative. Got: n = ",n)))
    
    g = TriangulatedPolygon(n)
    for i = 1:n-1
        add_edge!(g, i, i+1)
    end
    if n > 1
        add_edge!(g, n, 1)
    end
    if n <= 3
        return g
    end

    i = 2; j = n
    bo = true
    while i + 1 < j
        add_edge!(g, i, j)
        if bo
            i += 1
        else
            j -= 1
        end
        bo = !bo
    end
    return g
end


"""
    edges(g::TriangulatedPolygon) -> Vector{Edges}

Compute and return a list of all the edges in `g`.

Edges are not directed. It is however necessary to define a source and a target. 
For TriangulatedPolygon, the source will be the incident vertex with the smaller id.
"""
function edges(g::TriangulatedPolygon) :: Vector{Edge}
    E = collect(Edge(i,j) for i in 1:nv(g) for j in g.adjList[i])
    return filter!(e -> src(e) < dst(e), E)
end 

edgetype(g::TriangulatedPolygon) = SimpleEdge{Int}

"""
    has_edge(g::TriangulatedPolygon, e::Edge)
"""
has_edge(g::TriangulatedPolygon, e::Edge) = (dst(e) ∈ g.adjList[src(e)])
"""
    has_edge(g::TriangulatedPolygon, s, d)
"""
has_edge(g::TriangulatedPolygon, s, d) = (d ∈ g.adjList[s])

"""
    ne(g::TriangulatedPolygon)

Return the number of edges in `g`.
"""
has_vertex(g::TriangulatedPolygon, v) = (1 <= v <= g.n)

"""
    neighbors(g::TriangulatedPolygon, v)

Return a list of all the vertices in `g` that are adjacent to `v`.
"""
neighbors(g::TriangulatedPolygon, v) = g.adjList[v]
inneighbors(g::TriangulatedPolygon, v) = g.adjList[v]
outneighbors(g::TriangulatedPolygon, v) = g.adjList[v]

"""
    ne(g::TriangulatedPolygon)

Return the number of edges in `g`.
"""
ne(g::TriangulatedPolygon) = sum(size(g.adjList[i], 1) for i in 1:nv(g)) ÷ 2

"""
    nv(g::TriangulatedPolygon)

Return the number of vertices in `g`.
"""
nv(g::TriangulatedPolygon) = g.n

"""
    vertices(g::TriangulatedPolygon)

Return a list of all the vertices in `g`.
"""
vertices(g::TriangulatedPolygon) = collect(1:g.n)
is_directed(g::TriangulatedPolygon) = false
is_directed(::Type{TriangulatedPolygon}) = false

function add_edge!(g::TriangulatedPolygon, v, w) 
    if !has_edge(g, v, w)
        push!(g.adjList[v],w)
        push!(g.adjList[w],v)
    end
end

remove_edge!(g::TriangulatedPolygon, e::Edge) = remove_edge!(g,src(e),dst(e))
function remove_edge!(g::TriangulatedPolygon, src::Integer, dst::Integer)
    deleteat!(g.adjList[src], findfirst(x -> x==dst, g.adjList[src]))
    deleteat!(g.adjList[dst], findfirst(x -> x==src, g.adjList[dst]))
end


"""
    flip(g::TriangulatedPolygon, e::Edge)

Return the triangulated polygon obtained by flipping the edge `e` in `g`.
"""
flip(g::TriangulatedPolygon, e::Edge) = flip!(deepcopy(g), e) 

"""
    flip!(g::TriangulatedPolygon, src::Integer, dst::Integer)

Flip the the edge incident to `src` and `dst` in `g`.
"""
flip!(g::TriangulatedPolygon, src::Integer, dst::Integer) = flip!(g, Edge(src, dst))
"""
    flip!(g::TriangulatedPolygon, e::Edge)

Flip `e` in `g`.
"""
function flip!(g::TriangulatedPolygon, e::Edge)
    neigh1 = outneighbors(g, src(e))
    neigh2 = outneighbors(g, dst(e))
    S = intersect(neigh1, neigh2)
    u,v = S
    remove_edge!(g, e)
    add_edge!(g, u, v)
    return g
end

"""
    is_flippable(g::TriangulatedPolygon, e::Edge) -> Bool
"""
is_flippable(g::TriangulatedPolygon, e::Edge) = is_flippable(g, e.src, e.dst)
"""
    is_flippable(g::TriangulatedPolygon, src::Integer, dst::Integer) -> Bool

Return whether or not the edge can be flipped.

Note that for a triangulation of a convex polygon inner edges are always flippable, while outer edges cannot be flipped.    
"""
function is_flippable(g::TriangulatedPolygon, src::Integer, dst::Integer) :: Bool
    return length(intersect(outneighbors(g, src), outneighbors(g, dst))) >= 2
end


"""
    degrees(g::TriangulatedPolygon) -> Vector{Int}

Return a list of the degrees of every single vertex in `g`
"""
function degrees(g::TriangulatedPolygon) :: Vector{Int}
    return [length(g.adjList[i]) for i in 1:g.n]
end
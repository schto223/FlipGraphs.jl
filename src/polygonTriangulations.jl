"""
    struct TriangulatedPolygon <: AbstractGraph{Int32} 

A structure representing a triangulation of a convex polygon.
"""
struct TriangulatedPolygon{T<:Integer} <: AbstractGraph{T}
    n :: T
    adjList ::Vector{Vector{T}}

    function TriangulatedPolygon{T}(n::Integer) where T<:Integer
        new(n, Vector{Vector{T}}([[] for i in 1:n]))
    end
end

function Base.show(io::IO, mime::MIME"text/plain", g::TriangulatedPolygon)
    if g.n > 20
        println(io, string("TriangulatedPolygon with ", g.n, " vertices"))
    else
        println(io, string("TriangulatedPolygon with ", g.n, " vertices, and adjacency list:"))
        for i in eachindex(g.adjList)
            println(io, " $i $(i<10 ? " " : "")→ ["*join(g.adjList[i],", ")*"]")
        end
    end
end


"""
    triangulated_polygon(n::Integer) :: TriangulatedPolygon

Create a triangulated convex `n`-gon. 

Vertices are named from 1 to `n` in an anticlockwise manner.
The inside is triangulated in a zig-zag pattern.
"""
triangulated_polygon(n::Integer) :: TriangulatedPolygon = triangulated_polygon(typeof(n), n)
function triangulated_polygon(T::Type{<:Integer}, n::Integer) :: TriangulatedPolygon{T}
    n>=0 || throw(ArgumentError(string("The number of vertices(n) cannot be negative. Got: n = ",n)))
    
    #construct the polygon
    g = TriangulatedPolygon{T}(n)
    for i in 1:n-1
        add_edge!(g, i, i+1)
    end
    if n > 1
        add_edge!(g, n, 1)
    end
    if n <= 3
        return g
    end

    #add the inside edges in a zig-zag pattern
    i = 2; j = n
    bo = true
    while i + 1 < j
        add_edge!(g, i, j)
        bo ? i += 1 : j -= 1
        bo = !bo
    end
    return g
end


"""
    edges(g::TriangulatedPolygon) :: Vector{SimpleEdge{Int32}}

Compute and return a list of all the edges in `g`.

Edges are not directed. It is, however, necessary for computations to define a source and a target. 
For `TriangulatedPolygon`, the source will be the incident vertex with the smaller id.
"""
function edges(g::TriangulatedPolygon{T}) where T<:Integer
    return collect(SimpleEdge{T}, Edge(T(i),j) for i in 1:nv(g) for j in g.adjList[i] if i<j)
end 

edgetype(g::TriangulatedPolygon) = SimpleEdge{eltype(g)}

"""
    has_edge(g::TriangulatedPolygon, e::Edge)
"""
has_edge(g::TriangulatedPolygon, e::Edge) :: Bool = (dst(e) ∈ g.adjList[src(e)]) 
"""
    has_edge(g::TriangulatedPolygon, s::Integer, d::Integer)
    
Return `true` if `g` has an edge going from vertex `s` to vertex `d`. 
"""
has_edge(g::TriangulatedPolygon, s, d) = (d ∈ g.adjList[s])

"""
    has_vertex(g::TriangulatedPolygon, v::Integer)

Return `true` if `v` is a vertex in `g`. 
"""
has_vertex(g::TriangulatedPolygon, v) = (1 <= v <= g.n)

"""
    neighbors(g::TriangulatedPolygon, v::Integer) :: Vector{Int32}

Return the list of all the vertices in `g` that are adjacent to `v`.
"""
neighbors(g::TriangulatedPolygon, v::Integer) = g.adjList[v] 
inneighbors(g::TriangulatedPolygon, v) = g.adjList[v]
outneighbors(g::TriangulatedPolygon, v) = g.adjList[v]

"""
    ne(g::TriangulatedPolygon) :: Int

Return the number of edges in the triangulated convex polygon `g`.
"""
ne(g::TriangulatedPolygon) = sum(size(g.adjList[i], 1) for i in 1:nv(g)) ÷ 2

"""
    nv(g::TriangulatedPolygon) :: Int

Return the number of vertices/points in the triangulated convex polygon `g`.
"""
nv(g::TriangulatedPolygon) = g.n

"""
    vertices(g::TriangulatedPolygon) :: Vector{Int}

Create a list of all the vertices in the triangulated convex polygon `g`.
"""
vertices(g::TriangulatedPolygon) :: Vector{eltype(g)} = collect(eltype(g), 1:g.n)

is_directed(g::TriangulatedPolygon) = false
is_directed(::Type{TriangulatedPolygon}) = false

function add_edge!(g::TriangulatedPolygon, v, w)
    if !has_edge(g, v, w)
        push!(g.adjList[v],w)
        push!(g.adjList[w],v)
    end
end

remove_edge!(g::TriangulatedPolygon, e::Edge) = remove_edge!(g, src(e), dst(e))
function remove_edge!(g::TriangulatedPolygon, src::Integer, dst::Integer)
    deleteat!(g.adjList[src], findfirst(x -> x==dst, g.adjList[src]))
    deleteat!(g.adjList[dst], findfirst(x -> x==src, g.adjList[dst]))
end

"""
    flip(g::TriangulatedPolygon, src::Integer, dst::Integer) :: TriangulatedPolygon

Return the `TriangulatedPolygon` obtained from `g` by flipping the the edge incident to `src` and `dst`.
"""
flip(g::TriangulatedPolygon, src::Integer, dst::Integer) ::TriangulatedPolygon = flip!(deepcopy(g), src, dst)

"""
    flip(g::TriangulatedPolygon, e::Edge) :: TriangulatedPolygon

Return the triangulated polygon obtained by flipping the edge `e` in `g`.
"""
flip(g::TriangulatedPolygon, e::Edge) :: TriangulatedPolygon = flip(g, src(e), dst(e)) 

"""
    flip!(g::TriangulatedPolygon, e::Edge) :: TriangulatedPolygon

Flip the edge `e` in the triangulated convex polygon `g`.
"""
flip!(g::TriangulatedPolygon, e::Edge) :: TriangulatedPolygon = flip!(g, src(e), dst(e))

"""
    flip!(g::TriangulatedPolygon, src::Integer, dst::Integer) :: TriangulatedPolygon

Flip the the edge incident to `src` and `dst` in the triangulated convex polygon `g`.
"""
function flip!(g::TriangulatedPolygon, src::Integer, dst::Integer) :: TriangulatedPolygon
    u = 0; v = 0
    for j in g.adjList[src]
        if j in g.adjList[dst]
            if u == 0
                u = j
            else
                v = j
            end
        end
    end
    remove_edge!(g, src, dst)
    add_edge!(g, u, v)
    return g
end

"""
    flip_get_edge!(g::TriangulatedPolygon, src::Integer, dst::Integer) :: Tuple{Int32, Int32}

Flip the edge incident to the vertices `src` and `dst` and return the new endpoints of the flipped edge.
"""
function flip_get_edge!(g::TriangulatedPolygon{T}, src::Integer, dst::Integer) :: Tuple{T, T} where T<:Integer
    u = 0; v = 0
    for j in g.adjList[src]
        if j in g.adjList[dst]
            if u == 0
                u = j
            else
                v = j
            end
        end
    end
    remove_edge!(g, src, dst)
    add_edge!(g, u, v)
    return u,v
end

"""
    is_flippable(g::TriangulatedPolygon, e::Edge) :: Bool
"""
is_flippable(g::TriangulatedPolygon, e::Edge) ::Bool = is_flippable(g, e.src, e.dst)

"""
    is_flippable(g::TriangulatedPolygon, src::Integer, dst::Integer) :: Bool

Return whether or not the edge can be flipped.

Note that for a triangulation of a convex polygon, the inner edges are always flippable, 
    while the outer edges cannot be flipped.    
"""
function is_flippable(g::TriangulatedPolygon, src::Integer, dst::Integer) :: Bool
    return count(k in g.adjList[src] for k in g.adjList[dst]) ==2
end


"""
    degrees(g::TriangulatedPolygon) :: Vector{Int32}

Compute a list of the degrees of every single vertex in `g`.
"""
function degrees(g::TriangulatedPolygon{T}) :: Vector{T} where T<:Integer
    return [T(length(g.adjList[i])) for i in 1:g.n]
end


"""
    adjacency_matrix(g::TriangulatedPolygon) :: Matrix{Int32}

Compute the adjacency matrix for the triangulated graph `g`. 
"""
function adjacency_matrix(g::TriangulatedPolygon{T}) :: Matrix{T} where T<:Integer
    return adjacency_matrix(g.adjList)
end
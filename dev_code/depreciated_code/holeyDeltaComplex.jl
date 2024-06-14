
"""
    struct Crossing

A structure representing an edge passing through a hole.
"""
mutable struct Crossing
    hole_id::Int
    edge_id::Int
    going_in::Bool #true if the edge goes into the hole, false if it comes out
    
    key_holeposition::Int #the key to get the Crossing in the Hole \\TODO not sure if needed. remove if not used
    
    previous :: Crossing #previous crossing from the view of the hole
    next :: Crossing #next crossing from the view of the hole

    function Crossing(hole_id::Integer, edge_id::Integer, going_in::Bool)
        new(hole_id, edge_id, going_in, 0)
    end
end


"""
    struct Hole

A structure representing a **hole/handle** of an orientable surface.

It keeps track of the edges going in and out of the hole.\\
Although edges are not directed, we give them an arbitrary direction in order to distinct them.\\
The `Hole` may be imagined as a line going anticlockwise around the inner surface of the hole.\\
Edges, crossing this line from right to left are considered to be going into the hole, while edges crossing from left to right are coming out of the hole.
"""
mutable struct Hole
    id :: Int
    crossings :: Dict{Int, Crossing}
    last_key :: Int

    function Hole(id::Integer)
        new(id::Integer, Dict{Int, Crossing}(), 0)
    end
end

"""
    struct HoleyDeltaComplex

An extension of `DeltaComplex`, keeping track of which edges pass through which holes.

Many operations become slower compared to a `DeltaComplex`. 
However this structure is needed, if one wants to compare two different DeltaComplex'.
"""
struct HoleyDeltaComplex
    D::DeltaComplex
    holes::Vector{Hole}
    edge_crossings::Vector{Vector{Crossing}}
end

function Base.show(io::IO, mime::MIME"text/plain", C::Crossing)
    print(io, string("Crossing: (hole ",C.hole_id,")--<--",C.edge_id))
    if C.going_in
        print(io, "⤈") 
    else
        print(io, "⤉")
    end
    print(io, string("--<--"))
end

function Base.show(io::IO, C::Crossing)
    print(io, string(
        "Crossing(hole: ", C.hole_id,", key: ", C.key_holeposition, ", edge: ", C.edge_id,", going_in: ", C.going_in,
        ", previous: ", C.previous.key_holeposition,", next: ",C.next.key_holeposition,")"
        ))
end

function Base.show(io::IO, mime::MIME"text/plain", H::Hole)
    print(io, string("Hole ", H.id," : --<--"))
    C_first = H.crossings[minimum(keys(H.crossings))]
    if C_first.going_in
        print(io, string(C_first.edge_id, "⤈-"))
    else
        print(io, string(C_first.edge_id, "⤉-"))
    end    
    C = C_first.previous
    while C != C_first
        if C.going_in
            print(io, string(C.edge_id, "⤈-"))
        else
            print(io, string(C.edge_id, "⤉-"))
        end
        C = C.previous
    end
    print(io, string("-<--"))
end

function Base.show(io::IO, mime::MIME"text/plain", HD::HoleyDeltaComplex)
    print(io, string("HoleyDeltaComplex on orientable surface of genus ", genus(HD), " with: " ))
    print(io, string(np(HD), " Point", (np(HD)==1 ? "" : "s"),"; " ))
    println(io, string(nv(HD), " TriFaces","; " ))    
    print(io, string(ne(HD), " DualEdges","; " ))
    println(io, string(genus(HD), " Hole",(length(HD.holes)==1 ? "" : "s"),":"))
    for H in HD.holes
        show(io, mime, H); println(io)
    end
end


function insert_after!(H :: Hole, C_previous::Crossing, C::Crossing)
    H.last_key += 1
    H.crossings[H.last_key] = C
    C.key_holeposition = H.last_key
    C.next = C_previous.next
    C_previous.next.previous = C
    C_previous.next = C
    C.previous = C_previous
end

"""
    holeyDeltaComplex(g::Integer, num_points::Integer = 1) -> HoleyDeltaComplex

Create a **HoleyDeltaComplex** of genus `g`.
"""
function holey_delta_complex(g::Integer, num_points::Integer = 1) :: HoleyDeltaComplex
    function push_crossing!(H::Hole, C::Crossing)
        H.last_key += 1
        H.crossings[H.last_key] = C
        C.key_holeposition = H.last_key
        if H.last_key > 1
            H.crossings[H.last_key-1].next = C
            C.previous = H.crossings[H.last_key-1]
        end
    end

    if g == 0
        D = delta_complex(0,3)
    else  
        D = delta_complex(g)
    end

    holes = Vector{Hole}([Hole(i) for i in 1:genus(D)])
    edgeCrossings = Vector{Vector{Crossing}}([[] for i in 1:ne(D)])
    for i in 1:g
        if i <= g÷2 #holes on the upper side
            edges = [
                (get_edge(D, 6+8*(i-1), 3), false),
                (get_edge(D, 4+8*(i-1), 1), false),
                (get_edge(D, 4+8*(i-1), 3), false),
                (get_edge(D, 3+8*(i-1), 1), true),
            ]                
            if 2i!=g
                pushfirst!(edges, (get_edge(D, 6+8(i-1), 1), false))
            end
        else #holes on the lower side
            edges = [
                (get_edge(D, 1+8*(g-i), 3), true),
                (get_edge(D, 2+8*(g-i), 3), true),
                (get_edge(D, 2+8*(g-i), 1), true),
            ]                
            if 2*i-1 != g
                push!(edges, (get_edge(D, 4+8*(g-i), 3), true))
            end
            if 2*i-1 != g && i!=g
                push!(edges, (get_edge(D, 4+8*(g-i), 2), true))
            end
            if g==1 
                pop!(edges)
            end
        end
        for (e, going_in) in edges
            c = Crossing(i, e.id, going_in)
            push_crossing!(holes[i], c)
            push!(edgeCrossings[e.id], c)
        end
        H = holes[i]
        H.crossings[H.last_key].next = H.crossings[1]
        H.crossings[1].previous = H.crossings[H.last_key]
    end

    HD = HoleyDeltaComplex(D, holes, edgeCrossings)
    if g == 0
        num_points -= 2
    end
    for i in 2:num_points
        subdivide!(HD, i-1)
    end
    return HD
end

"""
    np(HD::HoleyDeltaComplex)

Return the number of points in the triangulation defined by `HD`.
"""
np(HD::HoleyDeltaComplex) = np(HD.D)

"""
    nv(HD::HoleyDeltaComplex)

Return the number of vertices(i.e. triangular faces) in `HD`.
"""
nv(HD::HoleyDeltaComplex) = nv(HD.D)

"""
    ne(HD::HoleyDeltaComplex)

Return the number of edges in the DeltaComplex `HD`.

This is equal to the number of edges in the triangulation itself.
"""
ne(HD::HoleyDeltaComplex) = ne(HD.D)

"""
    genus(HD::HoleyDeltaComplex) -> Int

Compute the genus of the HoleyDeltaComplex `HD`.
"""
genus(HD::HoleyDeltaComplex) ::Int = length(HD.holes)

"""
    euler_characteristic(HD::HoleyDeltaComplex) -> Int

Compute the euler characteristic of the HoleyDeltaComplex `HD`:\\
`` X = #vertices - #edges + #faces ``
"""
euler_characteristic(HD::HoleyDeltaComplex) :: Int = euler_characteristic(HD.D)

"""
    get_edge(HD::HoleyDeltaComplex, e::Integer) -> DualEdge
"""
get_edge(HD::HoleyDeltaComplex, e::Integer) = get_edge(HD.D, e)

"""
    get_edge(HD::HoleyDeltaComplex, t::Integer, side::Integer) -> DualEdge

Return the edge that forms the respective side in the given triangle `t`.
"""
get_edge(HD::HoleyDeltaComplex, t::Integer, side::Integer) = get_edge(HD.D, t, side)

"""
    get_vertex(HD::HoleyDeltaComplex, t::Integer) -> TriFace
"""
get_vertex(HD::HoleyDeltaComplex, t::Integer) = get_vertex(HD.D, t)

"""
    vertices(HD::HoleyDeltaComplex) -> Vector{TriFace}

Return the list of all vertices in `HD`.
"""
vertices(HD::HoleyDeltaComplex) = vertices(HD.D)

"""
    vertices(HD::HoleyDeltaComplex, d::DualEdge) -> Tuple{TriFace, TriFace}

Return both vertices adjacent to `d`.

The first TriFace is on the left of the edge, and the second one on the right of the edge.
"""
vertices(HD::HoleyDeltaComplex, d::DualEdge) :: Tuple{TriFace, TriFace} = vertices(HD.D, d)

"""
    edges(HD::HoleyDeltaComplex) -> Vector{DualEdge}

Return the list of all edges in `HD`. 
"""
edges(HD::HoleyDeltaComplex) = edges(HD.D)

"""
    edges(HD::HoleyDeltaComplex, t::Integer) -> Vector{DualEdge}

Return the list of all 3 edges that are incident to the `t`-th Triface in `HD`. 
"""
edges(HD::HoleyDeltaComplex, t::Integer) = edges(HD.D, t)

"""
    edges_id(HD::HoleyDeltaComplex, t::Integer) -> Tuple{Int, Int, Int}

Return the indices of all 3 edges that are incident to the `t`-th vertex. 
"""
edges_id(HD::HoleyDeltaComplex, t::Integer) = edges_id(HD.D, t)

"""
    points(HD::HoleyDeltaComplex, d::DualEdge) -> Tuple{Int, Int}

Return both endpoints of `d` in the triangulation.
"""
points(HD::HoleyDeltaComplex, d::DualEdge) = points(HD.D, d)


#is_orientable(HD::HoleyDeltaComplex) = is_orientable(HD.D)

"""
    point_degrees(HD::HoleyDeltaComplex) -> Vector{Int}

Return a vector containing the respective degree of each point in the triangulation.
"""
point_degrees(HD::HoleyDeltaComplex) ::Vector{Int} = point_degrees(HD.D)

"""
    diameter_triangulation(HD::HoleyDeltaComplex) -> Int

Compute the diameter of the triangulation defined by the HoleyDeltaComplex `HD`.\\

The **diameter** of a Graph is the greatest minimal distance between any 2 vertices.
"""
diameter_triangulation(HD::HoleyDeltaComplex) :: Int = diameter_triangulation(HD.D)

"""
    diameter_deltaComplex(HD::HoleyDeltaComplex)

Compute the diameter of the HoleyDeltaComplex `HD`.

The **diameter** of a Graph is the greatest minimal distance between any 2 vertices.

See also [`diameter_triangulation`](@ref)
"""
diameter_deltaComplex(HD::HoleyDeltaComplex) = diameter_deltaComplex(HD.D)

"""
    adjacency_matrix_deltaComplex(HD::HoleyDeltaComplex) :: Matrix{Int}

Compute the adjacency matrix of the delta complex `HD`.
"""
adjacency_matrix_deltaComplex(HD::HoleyDeltaComplex) = adjacency_matrix_deltaComplex(HD.D)

"""
    adjacency_matrix_triangulation(HD::HoleyDeltaComplex) :: Matrix{Int}

Compute the simple adjacency matrix of the triangulation defined by `HD`.

All entries are either 0 or 1.

See also [`multi_adjacency_matrix_triangulation`](@ref)
"""
adjacency_matrix_triangulation(HD::HoleyDeltaComplex) :: Matrix{Int} = adjacency_matrix_triangulation(HD.D)

"""
    multi_adjacency_matrix_triangulation(HD::HoleyDeltaComplex) :: Matrix{Int}

Compute the adjacency matrix of the multigraph of the triangulation defined by `HD`.

The ``i,j``-th entry not only notes if both points are connected, 
but also the number of edges that connect these two points.

See also [`adjacency_matrix_triangulation`](@ref)
"""
multi_adjacency_matrix_triangulation(HD::HoleyDeltaComplex) :: Matrix{Int} = multi_adjacency_matrix_triangulation(HD.D)

is_anticlockwise(HD::HoleyDeltaComplex, t::Integer, side::Integer) :: Bool = is_anticlockwise(HD.D, t, side)

"""
    triangle_edge(HD::HoleyDeltaComplex, d::DualEdge) -> Tuple{Int, Int}

Return the two points at either end of `d`.
"""
triangle_edge(HD::HoleyDeltaComplex, d::DualEdge) = triangle_edge(HD.D, d)

"""
    is_flippable(HD::HoleyDeltaComplex, e::Integer) -> Bool

Return true if the `e`-th edge can be flipped.

This is always the case if the edge does not connect a `TriFace` to itself.
"""
is_flippable(HD::HoleyDeltaComplex, e::Integer) :: Bool = is_flippable(HD.D, e)

"""
    flip(HD::HoleyDeltaComplex, d::DualEdge)

Return a copy of `HD` where the edge `d` has been flipped.
"""
function flip(HD::HoleyDeltaComplex, d::DualEdge; left::Bool = true) 
    return flip!(deepcopy(HD), d.id; left=left)
end

"""
    flip(HD::HoleyDeltaComplex, e::Integer)

Return a copy of `HD` where the `e`-th edge has been flipped.
"""
flip(HD::HoleyDeltaComplex, e::Integer; left::Bool = true) = flip!(deepcopy(HD), e, left=left)

"""
    flip!(HD::HoleyDeltaComplex, d::DualEdge)

Flip, if possible, the given edge in `HD`.
"""
flip!(HD::HoleyDeltaComplex, d::DualEdge; left::Bool = true) = flip!(HD, d.id, left=left)

"""
    flip!(HD::HoleyDeltaComplex, e::Integer)

Flip, if possible, the given edge in `HD`.
"""
function flip!(HD::HoleyDeltaComplex, e::Integer; left::Bool = true)
    d = get_edge(HD, e)  

    T1, T2 = vertices(HD.D, d)
    eA,eB,eC,eD = quadrilateral_edges(HD.D, d)
    ea,eb,ec,ed = eA.id, eB.id, eC.id, eD.id
    xy = is_anticlockwise(T1, right_side(get_side(d, 1)))
    yz = is_anticlockwise(T1, left_side(get_side(d, 1)))
    zq = is_anticlockwise(T2, right_side(get_side(d, 2)))
    qx = is_anticlockwise(T2, left_side(get_side(d, 2)))

    #remove edge crossings of e where the hole does no longer cross it after the edge is flipped
    i = 1
    while i <= length(edge_crossings(HD, e))
        c = edge_crossings(HD, e)[i]
        if c.going_in 
            if (c.next.edge_id == ed && c.next.going_in == qx) && (c.previous.edge_id == ea && c.previous.going_in == !xy) ||
                (c.next.edge_id == ec && c.next.going_in == zq) && (c.previous.edge_id == eb && c.previous.going_in == !yz)
                remove!(HD.holes[c.hole_id], c)
                deleteat!(HD.edge_crossings[e], i)
                i-=1
            end
        else
            if (c.previous.edge_id == ed && c.previous.going_in == !qx) && (c.next.edge_id == ea && c.next.going_in == xy) ||
                (c.previous.edge_id == ec && c.previous.going_in == !zq) && (c.next.edge_id == eb && c.next.going_in == yz)
                remove!(HD.holes[c.hole_id], c)
                deleteat!(HD.edge_crossings[e], i)
                i-=1
            end
        end
        i += 1
    end
    
    #reverse crossings along 
    reverse_edge_crossings = false
    for c in edge_crossings(HD, e)
        if c.next.edge_id == ea && c.next.going_in == xy || c.previous.edge_id == ea && c.previous.going_in == !xy
            reverse_edge_crossings = left
        else
            reverse_edge_crossings = !left
        end
    end
    if reverse_edge_crossings
        reverse!(HD.edge_crossings[e])
        for c in edge_crossings(HD, e)
            c.going_in = !c.going_in
        end
    end

    flip!(HD.D, e; left=left)
    
    # add crossing for edges that were parallel to e
    for c in (xy ? HD.edge_crossings[ea] : reverse(HD.edge_crossings[ea]))
        if c.going_in != xy #hole goes from ea to ...
            if c.next.edge_id == eb && (c.next.going_in == yz) #eb
                new_c = Crossing(c.hole_id, e, !left)
                insert_after!(HD.holes[c.hole_id], c, new_c)
                if left
                    push!(HD.edge_crossings[e], new_c)
                else 
                    pushfirst!(HD.edge_crossings[e], new_c)
                end
            end
        else #hole goes from ... to ea
            if c.previous.edge_id == eb && (c.previous.going_in == !yz)
                new_c = Crossing(c.hole_id, e, left)
                insert_after!(HD.holes[c.hole_id], c.previous, new_c)
                if left
                    push!(HD.edge_crossings[e], new_c)
                else 
                    pushfirst!(HD.edge_crossings[e], new_c)
                end
            end
        end
    end

    for c in (zq ? HD.edge_crossings[ec] : reverse(HD.edge_crossings[ec]))
        if c.going_in != zq #hole goes from ec to ...
            if c.next.edge_id == ed && (c.next.going_in == qx)
                new_c = Crossing(c.hole_id, e, left)
                insert_after!(HD.holes[c.hole_id], c, new_c)
                if !left
                    push!(HD.edge_crossings[e], new_c)
                else 
                    pushfirst!(HD.edge_crossings[e], new_c)
                end
            end
        else #hole goes from ... to ec
            if c.previous.edge_id == ed && (c.previous.going_in == !qx)
                new_c = Crossing(c.hole_id, e, !left)
                insert_after!(HD.holes[c.hole_id], c.previous, new_c)
                if !left
                    push!(HD.edge_crossings[e], new_c)
                else 
                    pushfirst!(HD.edge_crossings[e], new_c)
                end
            end
        end
    end
    return HD
end

"""
    edge_crossings(HD::HoleyDeltaComplex, e::Integer) -> Vector{Crossing}

Return an array of all the `Crossing`s involving the `e`-th edge.

The Crossings are ordered along the direction of `e`.
"""
function edge_crossings(HD::HoleyDeltaComplex, e::Integer) :: Vector{Crossing}
    return HD.edge_crossings[e]
end

"""
    hole_goes_left(HD::HoleyDeltaComplex, c::Crossing) -> Bool
    
Return `true` if the next edge crossing the hole is through the left side 
(relative to the entry point `c`) of the triangle.
"""
function hole_goes_left(HD::HoleyDeltaComplex, c::Crossing) :: Bool
    eNext_id = c.next.edge_id
    e = get_edge(HD, c.edge_id)
    side = c.going_in ? 2 : 1
    T = get_vertex(HD, e.triangles[side])
    left_e = get_edge(T, left_side(e.sides[side]))
    if left_e.id == eNext_id 
        if is_anticlockwise(T, left_side(e.sides[side])) == c.next.going_in
            return true
        end
    end 
    return false
end


function remove!(HD::HoleyDeltaComplex, c::Crossing)
    remove!(HD.holes[c.hole_id], c)
    filter!(cros -> c != cros, HD.edge_crossings[c.edge_id])
end

function remove!(H::Hole, c::Crossing)
    c.previous.next = c.next
    c.next.previous = c.previous
    delete!(H.crossings, c.key_holeposition)
end

"""
    remove_holeloops!(HD::HoleyDeltaComplex)

Remove all pairs of crossings where an edge enters a hole, and imeadiately comes out again.
"""
function remove_holeloops!(HD::HoleyDeltaComplex)
    foundLoop = false
    for H in HD.holes
        c_first = get_crossing(H)
        c = c_first
        c_last = c.previous
        i = length(H.crossings)
        while i > -1
            i -= 1 
            if c_last.edge_id == c.edge_id && c_last.going_in != c.going_in && (-1 <= indexin(c, HD.edge_crossings[c.edge_id])[1] - indexin(c_last, HD.edge_crossings[c.edge_id])[1] <= 1)
                edge_crossings = HD.edge_crossings[c.edge_id]
                delete!(edge_crossings, indexin(c, edge_crossings)[1])
                delete!(edge_crossings, indexin(c_last, edge_crossings)[1])
                remove!(HD, c)
                remove!(HD, c_last)
                c_last = c_last.previous
                c = c.next
                foundLoop = true
            else
                c_last = c
                c = c.next
            end
        end
    end
    if foundLoop
        remove_holeloops!(HD)
    end
    return HD
end

"""
    get_crossing(H::Hole) -> Crossing

Return any `Crossing` of the hole `H`.

Usually the `Crossing` which is returned, was the last one added to the `Hole`.
"""
function get_crossing(H::Hole) :: Crossing
    i = H.last_key
    while !(i in keys(H.crossings))
        i -= 1
    end
    return H.crossings[i]
end

"""
    subdivide!(HD ::HoleyDeltaComplex, t ::Integer)

Add a point to the inside of the `t`-th `TriFace` and connect it to each corner.
"""
function subdivide!(HD ::HoleyDeltaComplex, t ::Integer)
    e1,e2,e3 = edges_id(HD, t)  
    e4,e5,e6 = ne(HD)+1, ne(HD)+2, ne(HD)+3 #new edges that will be added in the inside of t

    t1,t2,t3 = t, nv(HD)+1, nv(HD)+2

    subdivide!(HD.D, t)
    
    push!(HD.edge_crossings,[],[],[])

    # holes passing through e1
    for c in (is_anticlockwise(HD, t1, 1) ? HD.edge_crossings[e1] : reverse(HD.edge_crossings[e1]))
        if !c.going_in == is_anticlockwise(HD, t1, 1) #the hole-line goes into the triangle
            if c.next.edge_id == e3
                c_new = Crossing(c.hole_id, e4, false)
                insert_after!(HD.holes[c.hole_id], c, c_new)
                push!(HD.edge_crossings[e4], c_new)
            else #c.next.edge_id == e2
                c_new = Crossing(c.hole_id, e5, true)
                insert_after!(HD.holes[c.hole_id], c, c_new)
                pushfirst!(HD.edge_crossings[e5], c_new)
            end
        else #c.going_in == is_anticlockwise(HD, T1.id, 1)
            if c.previous.edge_id == e3
                c_new = Crossing(c.hole_id, e4, true)
                insert_after!(HD.holes[c.hole_id], c.previous, c_new)
                push!(HD.edge_crossings[e4], c_new)
            else #c.previous.edge_id == e2
                c_new = Crossing(c.hole_id, e5, false)
                insert_after!(HD.holes[c.hole_id], c.previous, c_new)
                pushfirst!(HD.edge_crossings[e5], c_new)
            end
        end
    end

    # holes going from e2 to e3
    for c in (is_anticlockwise(HD, t2, 1) ? HD.edge_crossings[e2] : reverse(HD.edge_crossings[e2]))
        if !c.going_in == is_anticlockwise(HD, t2, 1) #the hole-line goes into the triangle
            if c.next.edge_id == e3
                c_new = Crossing(c.hole_id, e6, true)
                insert_after!(HD.holes[c.hole_id], c, c_new)
                pushfirst!(HD.edge_crossings[e6], c_new)
            end
        else #c.going_in == is_anticlockwise(HD, T1.id, 1)
            if c.previous.edge_id == e3
                c_new = Crossing(c.hole_id, e6, false)
                insert_after!(HD.holes[c.hole_id], c.previous, c_new)
                pushfirst!(HD.edge_crossings[e6], c_new)
            end
        end
    end

    return HD
end

"""
    num_crossings(H::Hole)

Return the number of `Crossing` that occur at `H`.
"""
function num_crossings(H::Hole) :: Integer
    c_first = get_crossing(H)
    c = c_first.next
    n = 1
    while c!=c_first
        c = c.next
        n += 1
    end
    return n
end

"""
    rename_points!(HD::HoleyDeltaComplex, p::Vector{<:Integer})

Rename all the points in `D` according to the permutation `p`.
"""
function rename_points!(HD::HoleyDeltaComplex, p::Vector{<:Integer})
    rename_points!(HD.D, p)
    return HD
end

"""
    rename_vertices!(HD::HoleyDeltaComplex, p::Vector{<:Integer})

Rename every vertex(TriFace) in `D`, according to the permutation `p`.

TriFace 1 => TriFace p[1]
"""
function rename_vertices!(HD::HoleyDeltaComplex, p::Vector{<:Integer})
    rename_vertices!(HD.D, p)
    return HD
end

"""
    rename_edges!(HD::HoleyDeltaComplex, p::Vector{<:Integer})
    
Rename every vertex(TriFace) in `HD`, according to the permutation `p`.

TriFace 1 => TriFace p[1]
"""
function rename_edges!(HD::HoleyDeltaComplex, p::Vector{<:Integer})
    rename_edges!(HD.D, p)
    for H in HD.holes
        foreach(c -> c.edge_id = p[c.edge_id], values(H.crossings))
    end
    HD.edge_crossings .= HD.edge_crossings[invert_permutation(p)]
    return HD
end


"""
    relative_point_degrees(HD::HoleyDeltaComplex, U::Vector{<:Integer}, V::Vector{<:Integer})

Return a vector containing the degree to `V` for each point in `U`.
"""
function relative_point_degrees(HD::HoleyDeltaComplex, U::Vector{<:Integer}, V::Vector{<:Integer})
    return relative_degrees(adjacency_matrix_triangulation(HD.D), U, V)
end


"""
    random_flips!(HD::HoleyDeltaComplex, n::Integer)

Randomly pick an edge, and flip it if possible. Repeat this `n` times.
"""
function random_flips!(HD::HoleyDeltaComplex, n::Integer)
    foreach(e-> flip!(HD, e), rand(HD.D.E, n))
    return D
end
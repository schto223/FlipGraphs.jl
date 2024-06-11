using StaticArrays
import Base.reverse, Base.reverse!

"""
    struct DualEdge

Representation of an edge in a DeltaComplex (i.e. the dual graph of a triangulation).

A DualEdge connects two TriFaces through specific sides.
"""
mutable struct DualEdge
    id :: Integer
    "The id's of the triangles incident to the DualEdge"
    triangles :: MVector{2, Int}
    "The respective sides on which the two triangles touch"
    sides :: MVector{2,Int8}  #the indices of the sides of the triangles through which the edge passes.
    is_twisted :: Bool     #true if one has to reverse one of the triangular faces.

    function DualEdge(t1::Integer, side1::Integer, t2::Integer, side2::Integer, is_twisted::Bool = false)
        new(0, (t1, t2), (side1, side2), is_twisted)
    end
end

"""
    struct TriFace

A `TriFace` represents a triangle in the triangulization of a surface. 

The TriFace's form the vertices of a DeltaComplex.
Each TriFace is formed between 3 points and is connected to 3 Trifaces through DualEdges.
The points and neighboring triangles do not have to be unique.

The points and edges are stored in an anticlockwise order.\\
The first edge/side is between the first and second point.\\
The second edge/side is between the second and third point.\\
The third edge/side is between the third and first point.\\
"""
struct TriFace
    id :: Base.RefValue{Int}  # the unique index number of this face in the DeltaComplex
    points :: MVector{3, Int} #corners x,y,z
    edges :: Vector{DualEdge}  #edge xy, yz, zx

    TriFace(id::Int, x::Int, y::Int, z::Int) = new(Ref(id), [x,y,z], Vector{DualEdge}(undef,3) )#MVector{3, DualEdge}(undef))
    TriFace(id::Int, points :: Vector{<:Integer}, edges:: Vector{DualEdge}) = new(Ref(id), points, edges)
end

"""
    struct DeltaComplex
    
A Graph datastructure representing a triangulation of a surface.

The DeltaComplex may be thought of as the dual of a triangulation.

Vertices are triangular faces (TriFace). Every vertex has three edges (DualEdge) incident to it.
"""
struct DeltaComplex
    V :: Vector{TriFace}
    E :: Vector{DualEdge}    
    num_points :: Base.RefValue{Int}

    function DeltaComplex()
        new(TriFace[], DualEdge[], Ref(1))
    end
end

function Base.show(io::IO, mime::MIME"text/plain", d::DualEdge)
    print(io, string("DualEdge ", d.id, " : (Δ",d.triangles[1],")-(",d.sides[1],")---"))
    if d.is_twisted
        print(io, "↺")
    else
        print(io, "-")
    end
    print(io, string("---(",d.sides[2],")-(Δ",d.triangles[2],")"))
end

function Base.show(io::IO, mime::MIME"text/plain", T::TriFace)
    print(io, string("TriFace #",id(T), ": Points(", T.points[1]," ", T.points[2]," ", T.points[3],")"))
    print(io, string(" Neighbors(",get_neighbor(T,1), " ",get_neighbor(T,2), " ",get_neighbor(T,3), ")"))
end

function Base.show(io::IO, mime::MIME"text/plain", D::DeltaComplex)
    if is_orientable(D)
        println(io, string("DeltaComplex on orientable surface of genus ", genus(D), " with ", np(D), " point",(np(D)==1 ? "" : "s")))
    else
        println(io, string("DeltaComplex on non-orientable surface of demigenus ", demigenus(D), " with ", np(D), " point", (np(D)==1 ? "" : "s")))
    end
    println(io, string(nv(D), " TriFaces:"))
    for T in vertices(D)
        print(io," "); show(io, mime, T); println(io)
    end
    println(io, string(ne(D), " DualEdges:"))
    for e in D.E
        print(io," "); show(io, mime, e); println(io)
    end  
end

"""
    vertices_id(d::DualEdge) -> Tuple{Int, Int}

Return the indices of both vertices (`TriFace`s) adjacent to `d`.
"""
vertices_id(d::DualEdge) :: Tuple{Int, Int} = d.triangles[1], d.triangles[2]

"""
    get_vertex_id(d::DualEdge, side::Integer) -> Int
 
Return the index of the vertex (`TriFace``) to the left of `d`.
"""
get_vertex_id(d::DualEdge, side::Integer) :: Int = d.triangles[side]

"""
    sides(d::DualEdge) -> Tuple{Int8, Int8}

Return the respective sides through which `d` connects the TriFaces.
"""
sides(d::DualEdge) ::Tuple{Int8, Int8} = d.sides[1], d.sides[2]

"""
    get_side(d::DualEdge, side::Integer) -> Int8

Return the respective side that `d` forms on the TriFace on the given `side`.
"""
get_side(d::DualEdge, side::Integer) :: Int8 = d.sides[side]

"""
    is_twisted(d::DualEdge) -> Bool

Return whether or not the given edge is twisted. 

If `d` is twisted, then everything on the other side gets regarded as its mirror image. 
"""
is_twisted(d::DualEdge) ::Bool = d.is_twisted

"""
    id(d::DualEdge) -> Int

Return the index of `d` in its `DeltaComplex`.
"""
id(d::DualEdge) :: Int = d.id

"""
    is_similar(d1::DualEdge, d2::DualEdge) -> Bool

Return `true` if d1 and d2 have the same twist and are connected to the same triangles.

This is only the case if `d1` and `d2` are the same edge, or if they are incident to a point of degree 2.
"""
function is_similar(d1::DualEdge, d2::DualEdge) :: Bool
    return d1.is_twisted == d2.is_twisted && (all(d1.triangles.==d2.triangles) || all(d1.triangles.==reverse(d2.triangles)))       
end

function update_endpoint!(d::DualEdge, i::Integer, ti::Integer, sidei::Integer)
    d.triangles[i] = ti
    d.sides[i] = sidei
end

#updates the endpoint of the dual Edge d
function update_endpoint!(d::DualEdge, t_old::Integer, side_old::Integer, t_new::Integer, side_new::Integer)
    if (d.triangles[1]==t_old && d.sides[1]==side_old)
        update_endpoint!(d, 1, t_new, side_new)
    else
        update_endpoint!(d, 2, t_new, side_new)
    end
end

"""
    other_endpoint(d::DualEdge, t::Integer, side::Integer) -> Tuple{Int, Int8}

Return the index of the other TriFace and its respective side, that is incident to `d` 
"""
function other_endpoint(d::DualEdge, t::Integer, side::Integer) :: Tuple{Int, Int8}
    if d.triangles[1] == t && d.sides[1] == side
        return d.triangles[2], d.sides[2]
    elseif d.triangles[2] == t && d.sides[2] == side
        return d.triangles[1], d.sides[1]
    else
        throw(ArgumentError(string("(t, side) = (", t,", ",side,") is not an endpoint of e = ",d)))
    end
end


has_point(T::TriFace, x::Integer) = (x in T.points) :: Bool

"""
    points(T::TriFace) -> Tuple{Int, Int, Int}

Return a tuple of the three points, that form the corners of `T`.
"""
points(T::TriFace) = Tuple(T.points) :: Tuple{Int, Int, Int}

"""
    get_point(T::TriFace, corner::Integer) -> Int
"""
get_point(T::TriFace, corner::Integer) = T.points[corner] :: Int

set_edge!(T::TriFace, side::Integer, d::DualEdge) = (T.edges[side] = d)

"""
    get_edge(T::TriFace, side::Integer) -> DualEdge
"""
get_edge(T::TriFace, side::Integer) = T.edges[side] :: DualEdge

"""
    get_edge_id(T::TriFace, side::Integer) -> Integer

Return the index of the edge on the respective side in `T`.
"""
get_edge_id(T::TriFace, side::Integer) = T.edges[side].id :: Integer

"""
    edges(T::TriFace) -> Vector{DualEdge}

Return the list of all 3 edges that are incident to `T`. 
"""
edges(T::TriFace) = T.edges :: Vector{DualEdge}

"""
    edges(D::DeltaComplex, t::Integer) -> Vector{DualEdge}

Return the list of all 3 edges that are incident to the `t`-th Triface in `D`. 
"""
edges(D::DeltaComplex, t::Integer) = D.V[t].edges :: Vector{DualEdge}

"""
    edges_id(D::DeltaComplex, t::Integer) -> Tuple{Int, Int, Int}

Return the indices of all 3 edges that are incident to `T`. 
"""
edges_id(D::DeltaComplex, t::Integer) = (get_edge_id(D.V[t],1), get_edge_id(D.V[t],2), get_edge_id(D.V[t],3)) :: Tuple{Int, Int, Int}

"""
    edges_id(T::TriFace) -> Tuple{Int, Int, Int}

Return the indices of all 3 edges that are incident to `T`. 
"""
edges_id(T::TriFace) = (T.edges[1].id, T.edges[2].id, T.edges[3].id) :: Tuple{Int, Int, Int}


"""
    id(T::TriFace)

Return the index of `T` in its DeltaComplex. 
"""
id(T::TriFace) = T.id.x

"""
    get_neighbor(T::TriFace, side::Integer) -> Int

Return the index of the TriFace adjacent to `T` on the given `side`.
"""
function get_neighbor(T::TriFace, side::Integer) ::Int
    if T.edges[side].triangles[1] == id(T)
        return T.edges[side].triangles[2]
    else
        return T.edges[side].triangles[1]
    end 
end

"""
    triangle_edge(T::TriFace, side::Integer) -> Tuple{Int, Int}

Return the two points incident to the respective side in the triangle `T`.
"""
function triangle_edge(T::TriFace, side::Integer) :: Tuple{Int, Int}
    return (T.points[side], T.points[(side%3) + 1])
end

"""
    add_vertex!(D::DeltaComplex, v::TriFace)
"""
add_vertex!(D::DeltaComplex, v::TriFace) = push!(D.V, v)

"""
    get_vertex(D::DeltaComplex, t::Integer) -> TriFace
"""
get_vertex(D::DeltaComplex, t::Integer)::TriFace = D.V[t]

"""
    vertices(D::DeltaComplex) -> Vector{TriFace}

Return the list of all vertices in `D`.
"""
vertices(D::DeltaComplex)::Vector{TriFace} = D.V 

"""
    vertices(D::DeltaComplex, d::DualEdge) -> Tuple{TriFace, TriFace}

Return both vertices adjacent to `d`.

The first TriFace is on the left of the edge, and the second one on the right of the edge.
"""
vertices(D::DeltaComplex, d::DualEdge) :: Tuple{TriFace, TriFace} = (D.V[get_vertex_id(d,1)], D.V[get_vertex_id(d,2)])

add_edge!(D::DeltaComplex, d::DualEdge) = push!(D.E, d)

"""
    get_edge(D::DeltaComplex, e::Integer) -> DualEdge
"""
get_edge(D::DeltaComplex, e::Integer) :: DualEdge = D.E[e] 

"""
    get_edge(D::DeltaComplex, t::Integer, side::Integer) -> DualEdge

Return the edge that forms the respective side in the given triangle `t`.
"""
get_edge(D::DeltaComplex, t::Integer, side::Integer) = D.V[t].edges[side]

"""
    edges(D::DeltaComplex) -> Vector{DualEdge}

Return the list of all the `DualEdge`s in `D`.
"""
edges(D::DeltaComplex) :: Vector{DualEdge} = D.E

set_num_points!(D::DeltaComplex, num_points) = setindex!(D.num_points, num_points)

"""
    points(D::DeltaComplex, d::DualEdge) -> Tuple{Int, Int}

Return both endpoints of `d` in the triangulation.
"""
points(D::DeltaComplex, d::DualEdge) = triangle_edge(D.V[d.triangles[1]], d.sides[1])

left_side(side::Integer) = (side+1)%3+1
right_side(side::Integer) = side%3 + 1


"""
    np(D::DeltaComplex)

Return the number of points in the triangulation defined by `D`.
"""
np(D::DeltaComplex) = getindex(D.num_points) :: Int

"""
    nv(D::DeltaComplex)

Return the number of vertices(i.e. triangular faces) in `D`.
"""
nv(D::DeltaComplex) = length(D.V) :: Integer

"""
    ne(D::DeltaComplex)

Return the number of edges in the DeltaComplex `D`.

This is equal to the number of edges in the triangulation itself.
"""
ne(D::DeltaComplex) = length(D.E) :: Integer


"""
    euler_characteristic(D::DeltaComplex)

Compute the euler characteristic of the DeltaComplex `D`:\\
`` X = #vertices - #edges + #faces ``
"""
euler_characteristic(D::DeltaComplex) = np(D) - ne(D) + nv(D) :: Int

"""
    genus(D::DeltaComplex)

Compute the genus of the DeltaComplex `D` if it forms an orientable surface.
"""
function genus(D::DeltaComplex)     
    if !is_orientable(D)
        throw(ArgumentError("Cannot compute the genus of a non-orientable surface"))
    end
    return (2-euler_characteristic(D))÷2  #assumes D is closed(i.e. has no boundary)
end

"""
    demigenus(D::DeltaComplex)

Compute the demigenus of the DeltaComplex `D`.\\

The **demigenus** or **non-orientable genus** (``k``) of a connected non-orientable closed surface is defined via the euler characteristic ()``X``) :\\
``k = 2 - X``
"""
function demigenus(D::DeltaComplex)
    if is_orientable(D)
        throw(ArgumentError("Cannot compute the demigenus of a orientable surface"))
    end
    return (2 - euler_characteristic(D))
end

"""
    adjacency_matrix_deltacomplex(D::DeltaComplex) :: Matrix{<:Integer}

Compute the adjacency matrix of the delta complex `D`.
"""
function adjacency_matrix_deltacomplex(D::DeltaComplex) :: Matrix{Int}
    A = zeros(Int, nv(D), nv(D))
    foreach(e-> (A[e.triangles[1], e.triangles[2]] = 1; A[e.triangles[2], e.triangles[1]] = 1), edges(D))
    return A
end

"""
    adjacency_matrix_triangulation(D::DeltaComplex) :: Matrix{Int}

Compute the simple adjacency matrix of the triangulation defined by `D`.

All entries are either 0 or 1.

See also [`multi_adjacency_matrix_triangulation`](@ref)
"""
function adjacency_matrix_triangulation(D::DeltaComplex) :: Matrix{Int}
    A = zeros(Int, np(D), np(D))
    foreach(t-> (
        A[t.points[1], t.points[2]] = 1; A[t.points[2], t.points[1]] = 1;
        A[t.points[2], t.points[3]] = 1; A[t.points[3], t.points[2]] = 1;
        A[t.points[3], t.points[1]] = 1; A[t.points[1], t.points[3]] = 1;
        ), vertices(D))
    return A
end

"""
    multi_adjacency_matrix_triangulation(D::DeltaComplex) :: Matrix{Int}

Compute the adjacency matrix of the multigraph of the triangulation defined by `D`.

The i,j-th entry not only notes if both points are connected, 
but also the number of edges that connect these two points.

See also [`adjacency_matrix_triangulation`](@ref)
"""
function multi_adjacency_matrix_triangulation(D::DeltaComplex) :: Matrix{Int}
    A = zeros(Int, np(D), np(D))
    foreach(t-> (
        A[t.points[1], t.points[2]] += 1; A[t.points[2], t.points[1]] += 1;
        A[t.points[2], t.points[3]] += 1; A[t.points[3], t.points[2]] += 1;
        A[t.points[3], t.points[1]] += 1; A[t.points[1], t.points[3]] += 1;
        ), vertices(D))
    return A.÷2
end

"""
    diameter_deltacomplex(D::DeltaComplex)

Compute the diameter of the DeltaComplex `D`.\\

The **diameter** of a Graph is the greatest minimal distance between any 2 vertices.

See also [`diameter_triangulation`](@ref)
"""
function diameter_deltaComplex(D::DeltaComplex)  
    return diameter(adjacency_matrix_deltacomplex(D))
end

"""
    diameter(D::DeltaComplex)

Compute the diameter of the DeltaComplex `D`.\\

The **diameter** of a Graph is the greatest minimal distance between any 2 vertices.

See also [`diameter_triangulation`](@ref)
"""
function diameter(D::DeltaComplex)
    return diameter_deltaComplex(D)
end


"""
    diameter_triangulation(D::DeltaComplex)

Compute the diameter of the triangulation defined by the DeltaComplex `D`.

The **diameter** of a Graph is the greatest minimal distance between any 2 vertices.

See also [`diameter_deltaComplex`](@ref)
"""
function diameter_triangulation(D::DeltaComplex)  
    return diameter(adjacency_matrix_triangulation(D))
end


"""
    glue_faces_along_edge!(D::DeltaComplex, t1:: TriFace, e1::Int, t2:: TriFace, e2::Int, twist::Bool=false)
    glue_faces_along_edge!(D::DeltaComplex, T1:: TriFace, e1::Integer, T2:: TriFace, e2::Integer, twist::Bool=false) 

Glue two triangle faces together along their shared edge. Only works if it is along the same edge in our triangulation.
"""
function glue_faces_along_edge!(D::DeltaComplex, t1:: Integer, e1::Integer, t2:: Integer, e2::Integer, twist::Bool=false)  
    return glue_faces_along_edge!(D, get_vertex(D, t1), e1, get_vertex(D, t2), e2, twist)
end
function glue_faces_along_edge!(D::DeltaComplex, T1:: TriFace, e1::Integer, T2:: TriFace, e2::Integer, twist::Bool=false) 
    x1,y1 = triangle_edge(T1, e1)
    x2,y2 = triangle_edge(T2, e2)

    if !twist
        x2,y2 = y2,x2
    end

    if x1 != x2
        merge_points!(D, x1, x2) 
        if y1 == max(x1,x2) 
            y1 = min(x1,x2)
        end
        if y2 == max(x1,x2)
            y2 = min(x1,x2)
        end
    end
    if y1 != y2
        merge_points!(D, y1, y2)
    end

    d = DualEdge(id(T1), e1, id(T2), e2, twist)
    add_edge!(D, d)
    set_edge!(T1, e1, d)
    set_edge!(T2, e2, d)
end

"""
    merge_points!(D::DeltaComplex, x::Integer, y::Integer)

Rename all references of the points `x` and `y` by their minimum.
"""
function merge_points!(D::DeltaComplex, x::Integer, y::Integer)
    if x > y
        x, y = y, x
    end
    if x != y
        set_num_points!(D, np(D)-1)
    end
    return rename_point!(D, y, x)
end

"""
    rename_point!(D::DeltaComplex, x_old::Integer, x_new::Integer)
"""
function rename_point!(D::DeltaComplex, x_old::Integer, x_new::Integer)
    foreach(T -> replace!(T.points, x_old => x_new), D.V)
    return D
end


"""
    deltacomplex(genus , num_points = 1)

Create a triangulation of an orientable surface with `num_points` points on it. 

By default `num_points` is set to 1.
"""
function deltacomplex(genus :: Integer, num_points :: Integer = 1)
    genus >= 0 || throw(ArgumentError(string("Cannot create a surface with a negative genus. Got: genus = ",genus)))
    genus != 0 || num_points >= 3 || throw(ArgumentError(string("To triangulate a spehere one needs at least 3 points.  Got: num_points = ", num_points))) 
    num_points>=1 || throw(ArgumentError(string("To triangulate a surface one needs at least 1 point.  Got: num_points = ", num_points))) 
    if genus == 0
        s = [1,-1,2,-2]
        num_points -= 2
    else
        n = 4*genus
        s = [(-1)^div(k-1, 2) * (2*((k-1)÷4) + (k-1)%2 + 1) for k in 1:n]
    end
    D = deltacomplex(s)
    for i in 2:num_points
        subdivide!(D, i-1)
    end
    return D
end

"""
    deltacomplex_non_orientable(demigenus , num_points = demigenus+1)

Create a triangulation of a non-orientable surface with `num_points` points on it. 

By default `num_points` is set to 1.
"""
function deltacomplex_non_orientable(demigenus ::Integer, num_points:: Integer = 1)
    demigenus > 0 || throw(ArgumentError(string("Cannot create a surface with a negative genus. Got: demigenus = ",demigenus)))
    num_points >= 1 || throw(ArgumentError(string("To triangulate a surface one needs at least 1 point.  Got: num_points = ", num_points))) 
    demigenus != 1 || num_points > 1 || throw(ArgumentError(string("To triangulate a projective plane one needs at least 2 points. Got: num_points = ", num_points))) 
    
    if demigenus % 2 == 0 #demigenus/2 klein bottles
        n = 2*demigenus
        s = [(-1)^(4%((k-1)%4+1)) * (2*((k-1)÷4) + (k-1)%2 + 1) for k in 1:n]
    elseif demigenus == 1
        s = [1,2,1,2] 
        num_points -= 1
    else                 #projective plane glued to (demigenus-1)/2 tori
        n = 2*(demigenus-3)
        s = [(-1)^div(k-1, 2) * (2*((k-1)÷4) + (k-1)%2 + 1) for k in 1:n]
        push!(s, demigenus-2, demigenus-1, demigenus, -demigenus+2, -demigenus+1, demigenus)
    end

    D = deltacomplex(s)
    for i in 1:(num_points-1)
        subdivide!(D, i)
    end
    return D
end

"""
    deltacomplex(s :: Array{<:Integer,1})

Create a triangulation of an orientable surface with a single point, by gluing corresponding edges together.\\

`s` should be an array of nonzero integers representing the edges of a polygon in anticlockwise order.\\ 
The i-th edge is orientated anticlockwise if `s[i]`>0 and anticlockwise if `s[i]`<0.\\
If `s[i]` and `s[j]` have the same absolute value, they are glued together while respecting their orientation.\\

#Examples
The following results in the triangulation of a torus with one point:
```julia-rep
julia> deltacomplex([1,2,-1,-2])
```

The following results in the triangulation of a Klein bottle with one point:
```julia-rep
julia> deltacomplex([1,2,-1,-2])
```
"""
function deltacomplex(s :: Array{<:Integer,1})
    n = length(s)

    function get_TriFace_and_rel_edge(D::DeltaComplex, i::Integer) :: Tuple{TriFace, Int8}
        if i < n÷2
            return get_vertex(D, 2*i-1), 1
        elseif i == n÷2
            return get_vertex(D, n-2), 1
        elseif i < n
            return get_vertex(D, 2*(n-i)), 2
        else
            return get_vertex(D, 1), 3
        end
    end

    s_abs = abs.(s)
    s_abs_unique = unique(s_abs)
    for i in s_abs_unique
        if count(==(i), s_abs) != 2
            if count(==(i), s_abs) == 1
                throw(ArgumentError(string(i)*" only appears once in s. \n
                    Delta-Complexes with boundary haven't been implemented.\n 
                Make sure, that every scaffolding edge is assigned to exactly 1 other edge."))
            else
                throw(ArgumentError(string(i)*" appears more than twice in s. \n
                    Each edge can only be assigned to 1 other edge"))
            end
        end
    end

    D = deltacomplex_scaffold(size(s,1))
    #glue respective sides of scaffold together
    for a in s_abs_unique
        i,j = findall(σ -> σ==a , s_abs)
        Ti, side_i = get_TriFace_and_rel_edge(D, i)
        Tj, side_j = get_TriFace_and_rel_edge(D, j)
        glue_faces_along_edge!(D, id(Ti), side_i, id(Tj), side_j, s[i]==s[j])
    end

    #rename points s.t. there are no "holes" from 1 to num_points
    pointnames = Set{Int}()
    foreach(T-> foreach(p-> push!(pointnames,p), points(T)), D.V)
    pnames = sort(collect(pointnames))
    for i in eachindex(pnames)
        if i != pnames[i]
            rename_point!(D,pnames[i],i)
        end
    end

    #give edges their id
    for i in eachindex(D.E)
        D.E[i].id = i
    end
    return D
end


function deltacomplex_scaffold(num_vertices :: Integer)
    num_vertices%2 == 0 || throw(ArgumentError(
        string("num_vertices has to be a multiple of 2. Got: num_vertices=", num_vertices)))
    num_vertices>0 || throw(ArgumentError(
        string("num_vertices has to be strictly positive. Got: num_vertices=", num_vertices)))
    
    n = num_vertices
    D = DeltaComplex()
    set_num_points!(D, n)

    T1 = TriFace(1, 1, 2, n)
    T2 = TriFace(2, 2, n-1, n)
    add_vertex!(D, T1)
    add_vertex!(D, T2)
    glue_faces_along_edge!(D, T1, 2, T2, 3)

    for i in 1 : (num_vertices-4) ÷ 2
        T1 = TriFace(2*i+1, i+1, i+2, n-i)
        add_vertex!(D, T1)
        glue_faces_along_edge!(D, T2, 1, T1, 3, false)
        T2 = TriFace(2*i+2, i+2, n-i-1, n-i)
        add_vertex!(D, T2)
        glue_faces_along_edge!(D, T1, 2, T2, 3, false)
    end
    return D
end


"""
    subdivide!(D::DeltaComplex, t::Integer)

Add a point to the inside of the `t`-th TriFace and connect it to each corner.
"""
function subdivide!(D::DeltaComplex, t::Integer)#twist
    set_num_points!(D, np(D)+1)
    n = nv(D)
    m = ne(D) 
    p = np(D)
    T = get_vertex(D, t)    
    v1,v2,v3 = points(T)
    d1,d2,d3 = edges(T)

    T1 = TriFace(id(T), [v1, v2, p], Vector{DualEdge}(undef,3)) # [get_neighbor(T,1) , n+1, n+2])
    T2 = TriFace(n+1, [v2, v3, p], Vector{DualEdge}(undef,3)) #[get_neighbor(T,2) , n+2, T.id])
    T3 = TriFace(n+2, [v3, v1, p], Vector{DualEdge}(undef,3)) #[get_neighbor(T,3) , T.id, n+1])
    D.V[t] = T1
    push!(D.V, T2, T3)
    
    d1_inside = (get_vertex_id(d1, 1) == id(T) && get_side(d1, 1) == 1) ? 1 : 2
    d2_inside = (get_vertex_id(d2, 1) == id(T) && get_side(d2, 1) == 2) ? 1 : 2
    d3_inside = (get_vertex_id(d3, 1) == id(T) && get_side(d3, 1) == 3) ? 1 : 2
    d1.sides[d1_inside] = 1
    d2.sides[d2_inside] = 1
    d3.sides[d3_inside] = 1
    d1.triangles[d1_inside] = id(T1)
    d2.triangles[d2_inside] = id(T2)
    d3.triangles[d3_inside] = id(T3)

    d4 = DualEdge(id(T3), 2, id(T1), 3, false)
    d4.id = m + 1
    d5 = DualEdge(id(T1), 2, id(T2), 3, false)
    d5.id = m + 2
    d6 = DualEdge(id(T2), 2, id(T3), 3, false)
    d6.id = m + 3
    push!(D.E,d4,d5,d6)

    T1.edges.=[d1, d5, d4]
    T2.edges.=[d2, d6, d5]
    T3.edges.=[d3, d4, d6]
    return D
end

"""
    twist_edges!(D::DeltaComplex, t::Integer)
    twist_edges!(T::TriFace)

Twist or untwist all 3 edges in a TriFace, and reverse the side order.

This action gives an equivalent representation ot the same triangulation.
It is usefull in the case that you would like to untwist a certain edge.

# Examples
```julia-repl
julia> D = deltacomplex([1,2,-1,2]);
julia> T = get_vertex(D,1);
julia> edges(T)
3-element Array{DualEdge,1}:
 DualEdge 2 : (Δ1)-(1)-------(2)-(Δ2)
 DualEdge 1 : (Δ1)-(2)-------(3)-(Δ2)
 DualEdge 3 : (Δ2)-(1)---↺---(3)-(Δ1)
 julia> twist_edges!(T);
 julia> edges(T)
 3-element Array{DualEdge,1}:
 DualEdge 3 : (Δ2)-(1)-------(1)-(Δ1)
 DualEdge 1 : (Δ1)-(2)---↺---(3)-(Δ2)
 DualEdge 2 : (Δ1)-(3)---↺---(2)-(Δ2)
```
"""
twist_edges!(D::DeltaComplex, t::Integer) = twist_edges!(get_vertex(D,t))
function twist_edges!(T::TriFace)
    update_endpoint!(get_edge(T,1), id(T), 1, id(T), 3)
    update_endpoint!(get_edge(T,3), id(T), 3, id(T), 1)
    reverse!(T.points)
    reverse!(T.edges)
    foreach(d -> d.is_twisted = !d.is_twisted, edges(T))
    return T
end


"""
    is_flippable(D::DeltaComplex, e::Integer)

Return true if the given edge can be flipped.

This is always the case if the edge does not connect a `TriFace` to itself.
"""
function is_flippable(D::DeltaComplex, e::Integer) 
    return is_flippable(get_edge(D,e))
end

"""
    is_flippable(d::DualEdge)

Return true if the given edge is can be flipped.

This is always the case if the edge does not connect a `TriFace` to itself.
"""
function is_flippable(d::DualEdge)
    return d.triangles[1] != d.triangles[2]
end

"""
    flip(D::DeltaComplex, e::Integer)

Return the resulting DeltaComplex obtained by flipping the given edge in `D`.
"""
function flip(D::DeltaComplex, e::Integer; left::Bool = true) 
    D2 = deepcopy(D)
    flip!(D2, get_edge(D2, e), left=left)
    return D2
end

"""
    flip!(D::DeltaComplex, e::Integer)

Flip, if possible, the given edge in `D`.
"""
flip!(D::DeltaComplex, e::Integer; left::Bool = true) = flip!(D, get_edge(D, e), left=left)

"""
    flip!(D::DeltaComplex, d::DualEdge)

Flip, if possible, the given edge in `D`.
"""
function flip!(D::DeltaComplex, e::DualEdge; left::Bool = true)
    is_flippable(e) || return false
    T1 = get_vertex(D, e.triangles[1])
    T2 = get_vertex(D, e.triangles[2])
    
    if e.is_twisted
        twist_edges!(T2)
    end 
    
    #find edges a,b,c,d and points x,y,z,q
    if e.sides[1] == 1
        t,a,b = edges(T1)
        z,x,y = points(T1)
        rot1 = 0
    elseif e.sides[1] == 2
        b,t,a = edges(T1)
        y,z,x = points(T1)
        rot1 = 1
    elseif e.sides[1] == 3
        a,b,t = edges(T1)
        x,y,z = points(T1)
        rot1 = -1
    end

    if e.sides[2] == 1
        t,c,d = edges(T2)
        q = get_point(T2,3)
        rot2 = 0
    elseif e.sides[2] == 2
        d,t,c = edges(T2)
        q = get_point(T2,1)
        rot2 = 1
    elseif e.sides[2] == 3
        c,d,t = edges(T2)
        q = get_point(T2,2)
        rot2 = -1
    end

    #The flip itself
    if left
        if e.sides[1] == 1
            T1.edges .= e,b,c
            T1.points .= q,y,z
        elseif e.sides[1] == 2
            T1.edges .= c,e,b
            T1.points .= z,q,y
        elseif e.sides[1] == 3
            T1.edges .= b,c,e
            T1.points .= y,z,q
        end

        if e.sides[2] == 1
            T2.edges .= e,d,a
            T2.points .= y,q,x
        elseif e.sides[2] == 2
            T2.edges .= a,e,d
            T2.points .= x,y,q
        elseif e.sides[2] == 3
            T2.edges .= d,a,e
            T2.points .= q,x,y
        end

        update_endpoint!(a, id(T1), 2 + rot1, id(T2), (3+rot2-1)%3+1)
        update_endpoint!(b, id(T1), (3+rot1-1)%3 + 1, id(T1), 2+rot1)
        update_endpoint!(c, id(T2), 2 + rot2, id(T1), (3+rot1-1)%3+1)
        update_endpoint!(d, id(T2), (3+rot2-1)%3 + 1, id(T2), 2+rot2)
    else #!left
        if e.sides[2] == 1
            T2.edges .= e,b,c
            T2.points .= q,y,z
        elseif e.sides[2] == 2
            T2.edges .= c,e,b
            T2.points .= z,q,y
        elseif e.sides[2] == 3
            T2.edges .= b,c,e
            T2.points .= y,z,q
        end

        if e.sides[1] == 1
            T1.edges .= e,d,a
            T1.points .= y,q,x
        elseif e.sides[1] == 2
            T1.edges .= a,e,d
            T1.points .= x,y,q
        elseif e.sides[1] == 3
            T1.edges .= d,a,e
            T1.points .= q,x,y
        end

        update_endpoint!(a, id(T1), 2 + rot1, id(T1), (3+rot1-1)%3+1)
        update_endpoint!(b, id(T1), (3+rot1-1)%3 + 1, id(T2), 2+rot2)
        update_endpoint!(c, id(T2), 2 + rot2, id(T2), (3+rot2-1)%3+1)
        update_endpoint!(d, id(T2), (3+rot2-1)%3 + 1, id(T1), 2+rot1)
    end
    return D
end

"""
    quadrilateral_edges(D::DeltaComplex, diagonal::DualEdge) -> Tuple{DualEdge, DualEdge, DualEdge, DualEdge}

Return the 4 edges who form a quadrilateral with th egiven `diagonal`.
"""
function quadrilateral_edges(D::DeltaComplex, diagonal::DualEdge) ::Tuple{DualEdge, DualEdge, DualEdge, DualEdge}
    t1,t2 = vertices_id(diagonal)
    return (
        get_edge(get_vertex(D,t1), right_side(diagonal.sides[1])),#a
        get_edge(get_vertex(D,t1), left_side(diagonal.sides[1])), #b
        get_edge(get_vertex(D,t2), right_side(diagonal.sides[2])),#c
        get_edge(get_vertex(D,t2), left_side(diagonal.sides[2]))  #d
    )
end

"""
    is_orientable(D::DeltaComplex)

Ckeck if the surface on which `D` lies is orientable or not.
"""
function is_orientable(D::DeltaComplex)
    color = zeros(Int8, nv(D))
    v0 = 1
    color[v0] = 1
    queue = [v0]
    while !isempty(queue)
        v = pop!(queue)
        for i in 1:3
            if color[get_neighbor(D.V[v], i)] == 0
                if  get_edge(D.V[v],i).is_twisted
                    color[get_neighbor(D.V[v], i)] = -color[v]
                else
                    color[get_neighbor(D.V[v], i)] = color[v]
                end
                push!(queue, get_neighbor(D.V[v], i))
            elseif (color[get_neighbor(D.V[v], i)] == color[v]) == (get_edge(D.V[v],i).is_twisted) 
                return false
            end
        end
    end
    return true
end

"""
    random_flips!(D::DeltaComplex, n::Integer)

Randomly pick an edge, and flip it if possible. Repeat this `n` times.
"""
function random_flips!(D::DeltaComplex, n::Integer)
    foreach(e-> flip!(D, e), rand(D.E,n))
    return D
end

"""
    randomize!(D::DeltaComplex; kwargs...) -> Int

Randomly flip edges in `D` until `D` is sufficiently generic.

Return the number of attempted flips.

The measure by which we determin if `D` is sufficiently generic is through its diameter.
This Method repeatedly flips a certain number of times.
After each flip sequence the diameter is computed.
Once this was repeated a certain number of times, the variance of all these past diameter measurements gets computed.

In theory, the variance should diminish over time. However, as we are flipping randomly, it will never truly converge to 0.
A certain flutter in the variance is expected, this will at some point cause the variance to increase every so often.
The algorithm stops once the last measured variance is bigger than the past few measurements.

# Arguments
- `num_initial_flips::Integer=1000000` : the number of flips to do before even start taking measurements.
- `num_flips_per_step::Integer=10000` : the number of flips to do before computing the diameter each step.
- `variance_interval_size::Integer=10` : the number of diameters to store, before computing their variance. 
- `lookback_size::Integer=2` : how far back to compare the current variance to, in order to decide when to stop. 

# Examples
```julia-repl
julia> D = deltacomplex(30,30);
julia> randomize!(D, num_initial_flips=10000, num_flips_per_step = 1000, variance_interval_size=10, lookback_size = 5)
160000
```
"""
function randomize!(D::DeltaComplex; num_initial_flips::Integer = 10000000, num_flips_per_step::Integer = 10000, variance_interval_size::Integer = 10, lookback_size::Integer = 2) :: Int
    function variance(x::Vector{<:Integer})
        mean = sum(x)//length(x)
        y = x.-mean 
        return sum(y.*y)//length(x)
    end
    n_flips = num_initial_flips +  variance_interval_size*num_flips_per_step
    pastvariances = [typemax(Int)//1 for i in 1:lookback_size]
    diameters = zeros(Int, variance_interval_size)
    random_flips!(D, num_initial_flips)
    diameters[1] = diameter_triangulation(D)
    
    for i in 2:variance_interval_size
        random_flips!(D,num_flips_per_step)
        diameters[i] = diameter_triangulation(D)
    end
    cur_variance = variance(diameters)
    k = 1
    while !all(cur_variance .>= pastvariances)
        pastvariances[k] = cur_variance
        k = k % lookback_size + 1
        for i in 1:variance_interval_size
            random_flips!(D, num_flips_per_step)
            diameters[i] = diameter_triangulation(D)
        end
        n_flips += variance_interval_size*num_flips_per_step
        cur_variance = variance(diameters)
    end
    return n_flips
end

"""
    point_degrees(D::DeltaComplex) -> Vector{<:Integer}

Return a vector containing the respective degree of each point in the triangulation.
"""
function point_degrees(D::DeltaComplex) ::Vector{<:Integer}
    pd = zeros(Int, np(D))
    foreach(T-> foreach(p->pd[p]+=1 , T.points) ,D.V)
    return pd
end


"""
    rename_points!(D::DeltaComplex, p::Vector{<:Integer})

Rename all the points in `D` according to the permutation `p`.
"""
function rename_points!(D::DeltaComplex, p::Vector{<:Integer})
    length(p) == np(D) || throw(ArgumentError("The permutation `p` does not have the right length.\n Expected vector of length $(np(D)) got: $(length(p))"))
    foreach(T -> T.points .= p[T.points], D.V)
    return D
end

"""
    rename_vertices!(D::DeltaComplex, p::Vector{<:Integer})

Rename every vertex(TriFace) in `D`, according to the permutation `p`.

TriFace 1 => TriFace p[1]
"""
function rename_vertices!(D::DeltaComplex, p::Vector{<:Integer})
    D.V .= D.V[invert_permutation(p)]
    foreach(T -> T.id.x = p[id(T)], D.V)
    foreach(d -> d.triangles .= p[d.triangles], D.E)
    return D
end

"""
    rename_edges!(D::DeltaComplex, p::Vector{<:Integer})

Rename every vertex(TriFace) in `D`, according to the permutation `p`.

TriFace 1 => TriFace p[1]
"""
function rename_edges!(D::DeltaComplex, p::Vector{<:Integer})
    D.E .= D.E[invert_permutation(p)]
    foreach(d -> d.id = p[d.id], D.E)
    return D
end
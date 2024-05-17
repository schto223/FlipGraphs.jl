using StaticArrays
import Base.reverse, Base.reverse!

export DeltaComplex, createDeltaComplex
export ne,np,nv
export get_edge, get_point, get_vertex, vertices, edges, points, edges_id, get_edge_id
export flip!, is_flippable, is_orientable, random_flips!, randomize!, point_degrees, relative_point_degrees
export euler_characteristic, genus, demigenus, diameter_triangulation, diameter_deltaComplex
export adjacency_matrix_triangulation, adjacency_matrix_deltaComplex
export subdivide!, twist_edges!
export rename_edges!, rename_points!, rename_trifaces!

mutable struct DualEdge
    id::Integer
    triangles :: MVector{2,Int}  #edge in dual of a triangulation
    sides :: MVector{2,Int8}  #the indices of the sides of the triangles through which the edge passes.
    is_twisted :: Bool     #true if the one has to flip one of the triangular faces.

    function DualEdge(t1::Integer, side1::Integer, t2::Integer, side2::Integer, is_twisted::Bool = false)
        new(0, (t1, t2), (side1, side2), is_twisted)
    end
    function DualEdge(triangles :: MVector{2, <:Integer}, sides::MVector{2, <:Integer}, is_twisted::Bool = false)
        new(0, triangles, sides, is_twisted)
    end
end

struct TriFace
    id :: Int  # the unique index number of this face in the DeltaComplex
    points :: MVector{3, Int} #corners x,y,z
    edges :: Vector{DualEdge}  #edge xy, yz, zx

    TriFace(id::Int, x::Int, y::Int, z::Int) = new(id, [x,y,z], Vector{DualEdge}(undef,3) )#MVector{3, DualEdge}(undef))
    TriFace(id::Int, points :: Vector{<:Integer}, edges:: Vector{DualEdge}) = new(id,points, edges)
end

"""
    struct DeltaComplex
    
A Graph datastructure representing a triangulation of a surface. \\
Vertices are triangular faces. Every vertex has three edges incident to it.
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
    print(io, string("TriFace #",T.id, ": Points(", T.points[1]," ", T.points[2]," ", T.points[3],")"))
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

vertices(d::DualEdge) = d.triangles[1], d.triangles[2]
get_vertex(d::DualEdge, i::Integer) = d.triangles[i]
sides(d::DualEdge) = d.sides[1], d.sides[2]
get_side(d::DualEdge, i::Integer) = d.sides[i]
set_twisted!(d:: DualEdge, is_twisted::Bool) = (d.is_twisted = is_twisted)
is_twisted(d::DualEdge) = d.is_twisted
get_id(d::DualEdge) = d.id

reverse(d::DualEdge) = DualEdge(reverse(d.triangles), reverse(d.sides), d.is_twisted)
function reverse!(d::DualEdge) 
    reverse!(d.triangles)
    reverse!(d.sides)
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

function other_endpoint(e::DualEdge, t::Integer, side::Integer) :: Tuple{Int, Int8}
    if e.triangles[1] == t && e.sides[1] == side
        return e.triangles[2], e.sides[2]
    elseif e.triangles[2] == t && e.sides[2] == side
        return e.triangles[1], e.sides[1]
    else
        throw(ArgumentError(string("(t, side) = (", t,", ",side,") is not an endpoint of e = ",e)))
    end
end


has_point(T::TriFace, x::Integer) = (x in T.points) :: Bool
points(T::TriFace) = return Tuple(T.points) :: Tuple{Int, Int, Int}
get_point(T::TriFace, corner::Integer) = T.points[corner]
set_edge!(T::TriFace, side::Integer, edge::DualEdge) = (T.edges[side] = edge)
get_edge(T::TriFace, side::Integer) = T.edges[side] :: DualEdge
get_edge_id(T::TriFace, side::Integer) = T.edges[side].id :: Integer
edges(T::TriFace) = T.edges :: Vector{DualEdge}
edges(D::DeltaComplex, t::Integer) = D.V[t].edges :: Vector{DualEdge}
edges_id(D::DeltaComplex, t::Integer) = (get_edge_id(D.V[t],1), get_edge_id(D.V[t],2), get_edge_id(D.V[t],3)) :: Tuple{Int, Int, Int}

function get_neighbor(T::TriFace, side::Integer) ::Int
    if T.edges[side].triangles[1] == T.id
        return T.edges[side].triangles[2]
    else
        return T.edges[side].triangles[1]
    end 
end

function triangle_edge(T::TriFace, side::Integer)
    return (T.points[side], T.points[(side%3) + 1])
end

is_anticlockwise(D::DeltaComplex, t::Integer, side::Integer) :: Bool = is_anticlockwise(D.V[t], side)
function is_anticlockwise(T::TriFace, side::Integer) :: Bool
    d = get_edge(T, side)
    return (T.id == get_vertex(d, 1) && get_side(d,1) == side)
end


add_vertex!(D::DeltaComplex, v::TriFace) = push!(D.V, v)
get_vertex(D::DeltaComplex, i::Integer)::TriFace = D.V[i]
vertices(D::DeltaComplex)::Vector{TriFace} = D.V 
vertices(D::DeltaComplex, d::DualEdge) :: Tuple{TriFace, TriFace} = (D.V[get_vertex(d,1)], D.V[get_vertex(d,2)])
add_edge!(D::DeltaComplex, d::DualEdge) = push!(D.E, d)
get_edge(D::DeltaComplex, i::Integer) = D.E[i]
get_edge(D::DeltaComplex, t::Integer, side::Integer) = D.V[t].edges[side]
edges(D::DeltaComplex) = D.E
remove_edge!(D::DeltaComplex, d::DualEdge) = remove!(D.E, d)
set_num_points!(D::DeltaComplex, num_points) = setindex!(D.num_points, num_points)
points(D::DeltaComplex, d::DualEdge) = triangle_edge(D.V[d.triangles[1]], d.sides[1])

left_side(side::Integer) = (side+1)%3+1
right_side(side::Integer) = side%3 + 1

"""
    left_edge(D::DeltaComplex, t::Integer, side::Integer) :: (DualEdge, Int, Int8)

Return the edge to the left of the side `side` in the `t`-th triangle-vertex of `D`, \\
the other endpoint and \\
its side in the other endpoint. \\
"""
function left_edge(D::DeltaComplex, t::Integer, side::Integer) :: Tuple{DualEdge, Int, Int8}
    side_out = left_side(side)
    e = get_edge(D, t, side_out)
    t_new, side_in = other_endpoint(e, t, side_out)
    return e, t_new, side_in
end

function right_edge(D::DeltaComplex, t::Integer, side::Integer) :: Tuple{DualEdge, Int, Int8}
    side_out = right_side(side)
    e = get_edge(D, t, side_out)
    t_new, side_in = other_endpoint(e, t, side)
    return e, t_new, side_in
end

"""
    np(D::DeltaComplex)

Return the number of points in the triangulation defined by `D`
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
    adjacency_matrix_deltaComplex(D::DeltaComplex) :: Matrix{<:Integer}

Compute the adjacency matrix of the delta complex `D`.
"""
function adjacency_matrix_deltaComplex(D::DeltaComplex) :: Matrix{<:Integer}
    A = zeros(Int, nv(D), nv(D))
    foreach(e-> (A[e.triangles[1], e.triangles[2]] = 1; A[e.triangles[2], e.triangles[1]] = 1), edges(D))
    return A
end

"""
    adjacency_matrix_triangulation(D::DeltaComplex) :: Matrix{<:Integer}

Compute the adjacency matrix of the triangulation defined by `D`.
"""
function adjacency_matrix_triangulation(D::DeltaComplex) :: Matrix{<:Integer}
    A = zeros(Int, np(D), np(D))
    foreach(t-> (
        A[t.points[1], t.points[2]] = 1; A[t.points[2], t.points[1]] = 1;
        A[t.points[2], t.points[3]] = 1; A[t.points[3], t.points[2]] = 1;
        A[t.points[3], t.points[1]] = 1; A[t.points[1], t.points[3]] = 1;
        ), vertices(D))
    return A
end


"""
    diameter_deltaComplex(D::DeltaComplex)

Compute the diameter of the DeltaComplex `D`.\\

The **diameter** of a Graph is the greatest minimal distance between any 2 vertices.
"""
function diameter_deltaComplex(D::DeltaComplex)  
    A = adjacency_matrix_deltaComplex(D)
    return diameter(A)
end


"""
    diameter_triangulation(D::DeltaComplex)

Compute the diameter of the triangulation defined by the DeltaComplex `D`.\\

The **diameter** of a Graph is the greatest minimal distance between any 2 vertices.
"""
function diameter_triangulation(D::DeltaComplex)  
    A = adjacency_matrix_triangulation(D)
    return diameter(A)
end


"""
    glue_faces_along_edge!(D::DeltaComplex, t1:: TriFace, edge_index_1::Int, t2:: TriFace, edge_index_2::Int)

Glue two triangle faces together along their shared edge. Only works if it is along the same edge in our triangulation.
"""
function glue_faces_along_edge!(D::DeltaComplex, t1:: Integer, edge_index_1::Integer, t2:: Integer, edge_index_2::Integer, twist::Bool=false)  
    return glue_faces_along_edge!(D, get_vertex(D, t1), edge_index_1, get_vertex(D, t2), edge_index_2, twist)
end

function glue_faces_along_edge!(D::DeltaComplex, T1:: TriFace, edge_index_1::Integer, T2:: TriFace, edge_index_2::Integer, twist::Bool=false) 

    #set_neighbor!(T1, edge_index_1, t2)
    #set_neighbor!(T2, edge_index_2, t1)

    x1,y1 = triangle_edge(T1, edge_index_1)
    x2,y2 = triangle_edge(T2, edge_index_2)

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

    e = DualEdge(T1.id, edge_index_1, T2.id, edge_index_2, twist)
    add_edge!(D,e)
    set_edge!(T1, edge_index_1, e)
    set_edge!(T2, edge_index_2, e)
end


function merge_points!(D::DeltaComplex, x::Integer, y::Integer)
    if x > y
        x, y = y, x
    end
    if x!=y
        set_num_points!(D, np(D)-1)
    end
    return rename_point!(D, y, x)
end


function rename_point!(D::DeltaComplex, x_old::Integer, x_new::Integer)
    foreach(T -> replace!(T.points, x_old => x_new), D.V)
    return D
end


"""
    createDeltaComplex(genus [, num_points])

Create a triangulation of an orientable surface with `num_points` points on it. \\
By default `num_points` is set to 1.
"""
function createDeltaComplex(genus :: Integer, num_points :: Integer = 1)
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
    D = createDeltaComplex(s)
    for i in 2:num_points
        subdivide!(D, i-1)
    end
    return D
end

"""
    createDeltaComplex(s :: Array{<:Integer,1})

Create a triangulation of an orientable surface with a single point, by gluing corresponding edges together.\\

`s` should be an array of nonzero integers representing the edges of a polygon in anticlockwise order.\\ 
The i-th edge is orientated anticlockwise if `s[i]`>0 and anticlockwise if `s[i]`<0.\\
If `s[i]` and `s[j]` have the same absolute value, they are glued together while respecting their orientation.\\

#Examples
The following results in the triangulation of a torus with one point:
```julia-rep
julia> createDeltaComplex([1,2,-1,-2])
```

The following results in the triangulation of a Klein bottle with one point:
```julia-rep
julia> createDeltaComplex([1,2,-1,-2])
```
"""
function createDeltaComplex(s :: Array{<:Integer,1})
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

    D = createDeltaScaffold(size(s,1))
    #glue respective sides of scaffold together
    for a in s_abs_unique
        i,j = findall(σ -> σ==a , s_abs)
        Ti, side_i = get_TriFace_and_rel_edge(D, i)
        Tj, side_j = get_TriFace_and_rel_edge(D, j)
        glue_faces_along_edge!(D, Ti.id, side_i, Tj.id, side_j, s[i]==s[j])
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


function createDeltaScaffold(num_vertices :: Integer)
    num_vertices%2 == 0 || throw(ArgumentError(
        string("num_vertices has to be a multiple of 2. Got: num_vertices=",num_vertices)))
    num_vertices>0 || throw(ArgumentError(
        string("num_vertices has to be strictly positive. Got: num_vertices=",num_vertices)))
    
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

    T1 = TriFace(T.id, [v1, v2, p], Vector{DualEdge}(undef,3)) # [get_neighbor(T,1) , n+1, n+2])
    T2 = TriFace(n+1, [v2, v3, p], Vector{DualEdge}(undef,3)) #[get_neighbor(T,2) , n+2, T.id])
    T3 = TriFace(n+2, [v3, v1, p], Vector{DualEdge}(undef,3)) #[get_neighbor(T,3) , T.id, n+1])
    D.V[t] = T1
    push!(D.V, T2, T3)
    
    d1_inside = (get_vertex(d1, 1) == T.id && get_side(d1, 1) == 1) ? 1 : 2
    d2_inside = (get_vertex(d2, 1) == T.id && get_side(d2, 1) == 2) ? 1 : 2
    d3_inside = (get_vertex(d3, 1) == T.id && get_side(d3, 1) == 3) ? 1 : 2
    d1.sides[d1_inside] = 1
    d2.sides[d2_inside] = 1
    d3.sides[d3_inside] = 1
    d1.triangles[d1_inside] = T1.id
    d2.triangles[d2_inside] = T2.id
    d3.triangles[d3_inside] = T3.id

    d4 = DualEdge(T3.id, 2, T1.id, 3, false)
    d4.id = m+1
    d5 = DualEdge(T1.id, 2, T2.id, 3, false)
    d5.id = m+2
    d6 = DualEdge(T2.id, 2, T3.id, 3, false)
    d6.id = m+3
    push!(D.E,d4,d5,d6)

    T1.edges.=[d1,d5,d4]
    T2.edges.=[d2,d6,d5]
    T3.edges.=[d3,d4,d6]
    return D
end

twist_edges!(D::DeltaComplex, t::Integer) = twist_edges!(get_vertex(D,t))
function twist_edges!(T::TriFace)
    update_endpoint!(get_edge(T,1), T.id, 1, T.id, 3)
    update_endpoint!(get_edge(T,3), T.id, 3, T.id, 1)
    reverse!(T.points)
    reverse!(T.edges)
    foreach(d -> d.is_twisted = !d.is_twisted, edges(T))
    return T
end

"""
    is_flippable(D::DeltaComplex, e::Integer)
    is_flippable(d::DualEdge)

Return true if the given edge is can be flipped.\\ 
This is only the case if the edge does not connect a triangle face to itself.
"""
is_flippable(D::DeltaComplex, e::Integer) = is_flippable(get_edge(D,e))
function is_flippable(d::DualEdge)
    return d.triangles[1] != d.triangles[2]
end

"""
    flip!(D::DeltaComplex, e::Integer)
    flip!(D::DeltaComplex, e::DualEdge)

Flip the given edge in `D`
"""
flip!(D::DeltaComplex, e::Integer, left::Bool = true) = flip!(D, get_edge(D, e), left)
function flip!(D::DeltaComplex, e::DualEdge, left::Bool = true) :: Bool
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
    #if !left
    #    T1 = get_vertex(D, e.triangles[2])
    #    T2 = get_vertex(D, e.triangles[1])
    #end

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

        update_endpoint!(a, T1.id, 2 + rot1, T2.id, (3+rot2-1)%3+1)
        update_endpoint!(b, T1.id, (3+rot1-1)%3 + 1, T1.id, 2+rot1)
        update_endpoint!(c, T2.id, 2 + rot2, T1.id, (3+rot1-1)%3+1)
        update_endpoint!(d, T2.id, (3+rot2-1)%3 + 1, T2.id, 2+rot2)
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

        update_endpoint!(a, T1.id, 2 + rot1, T1.id, (3+rot1-1)%3+1)
        update_endpoint!(b, T1.id, (3+rot1-1)%3 + 1, T2.id, 2+rot2)
        update_endpoint!(c, T2.id, 2 + rot2, T2.id, (3+rot2-1)%3+1)
        update_endpoint!(d, T2.id, (3+rot2-1)%3 + 1, T1.id, 2+rot1)
    end
    return true
end

"""
    quadrilateral_edges(D::DeltaComplex, diagonal::DualEdge) ::Tuple{Integer, Integer, Integer, Integer}
"""
function quadrilateral_edges(D::DeltaComplex, diagonal::DualEdge) ::Tuple{DualEdge, DualEdge, DualEdge, DualEdge}
    t1,t2 = vertices(diagonal)
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
                    color[get_neighbor(D.V[v], i)] = - color[v]
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


function randomize!(D::DeltaComplex, n_flips_initial::Integer, n_flips_step::Integer, variance_interval_size::Integer)
    function variance(x::Vector{<:Integer})
        mean = sum(x)//length(x)
        y = x.-mean 
        return sum(y.*y)//length(x)
    end
    
    diameters = zeros(Int, variance_interval_size)
    random_flips!(D,n_flips_initial)
    diameters[1] = diameter_triangulation(D)
    
    for i in 2:variance_interval_size
        random_flips!(D,n_flips_step)
        diameters[i] = diameter_triangulation(D)
    end
    cur_variance = variance(diameters)
    last_variance = cur_variance + 1
    while cur_variance < last_variance
        last_variance = cur_variance
        for i in 1:variance_interval_size
            random_flips!(D,n_flips_step)
            diameters[i] = diameter_triangulation(D)
        end
        cur_variance = variance(diameters)
    end
    return D
end

"""
    point_degrees(D::DeltaComplex)

Return a vector containing the respective degree of each point in the triangulation.
"""
function point_degrees(D::DeltaComplex) ::Vector{<:Integer}
    pd = zeros(Int, np(D))
    foreach(T-> foreach(p->pd[p]+=1 , T.points) ,D.V)
    return pd
end

"""
    relative_point_degrees(D::DeltaComplex, U::Vector{<:Integer}, V::Vector{<:Integer})

Return a vector containing the degree to `V` for each point in `U`.
"""
function relative_point_degrees(D::DeltaComplex, U::Vector{<:Integer}, V::Vector{<:Integer})
    A = adjacency_matrix_triangulation(D)
    return [sum(A[u,V]) for u in U]
end

function rename_points!(D::DeltaComplex, p::Vector{<:Integer})
    foreach(T -> T.points = p[T.points], D.V)
    return D
end

function rename_trifaces!(D::DeltaComplex, p::Vector{<:Integer})
    D.V = D.V[invert_perm(p)]
    foreach(T -> T.id = p[T.id], D.V)
    foreach(d -> d.triangles = p[d.triangles], D.E)
    return D
end

function rename_edges!(D::DeltaComplex, p::Vector{<:Integer})
    D.E = D.E[invert_perm(p)]
    foreach(d -> d.id = p[d.id], D.E)
    return D
end
using StaticArrays
import Base.reverse

export DeltaComplex, createDeltaComplex,
get_num_edges, get_num_points, get_num_trifaces, ne,np,nv,
get_edge, get_points, get_vertex,
trifaces, flip!, is_flippable, is_orientable,
euler_characteristic, genus, 
adjacency_matrix

struct DualEdge
    triangles :: MVector{2,Int}  #edge in dual of a triangulation
    resp_sides :: MVector{2, Int8}  #the indices of the sides of the triangles through which the edge passes.
    is_twisted :: Bool     #true if the one has to flip one of the triangular faces.

    function DualEdge(t1::Integer, side1::Integer, t2::Integer, side2::Integer, is_twisted::Bool = false)
        new((t1, t2), (side1, side2), is_twisted)
    end
    function DualEdge(triangles :: MVector{2, <:Integer}, resp_sides::MVector{2, <:Integer}, is_twisted::Bool)
        new(triangles, resp_sides, is_twisted)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", d::DualEdge)
    print(io, string("DualEdge: (Δ",d.triangles[1],")-(",d.resp_sides[1],")---"))
    if d.is_twisted
        print(io, "↺")
    else
        print(io, "-")
    end
    print(io, string("---(",d.resp_sides[2],")-(Δ",d.triangles[2],")"))
end

function update_endpoints!(d::DualEdge, t1::Integer, side1::Integer, t2::Integer, side2::Integer)
    d.triangles .= t1,t2
    d.resp_sides .= side1,side2
end

function update_endpoint!(d::DualEdge, i::Integer, ti::Integer, sidei::Integer)
    d.triangles[i] = ti
    d.resp_sides[i] = sidei
end

function update_endpoint!(d::DualEdge, t_old::Integer, side_old::Integer, t_new::Integer, side_new::Integer)
    if (d.triangles[1]==t_old && d.resp_sides[1]==side_old)
        update_endpoint!(d, 1, t_new, side_new)
    else
        update_endpoint!(d, 2, t_new, side_new)
    end
end

reverse(d::DualEdge) = DualEdge(reverse(d.triangles), reverse(d.resp_sides), d.is_twisted)
is_twisted(d::DualEdge) = d.is_twisted

function get_other_endpoint(d::DualEdge, t::Integer) :: Int
    if d.triangles[1] == t
        return d.triangles[2]
    else
        return d.triangles[1]
    end
end

function set_twisted!(e:: DualEdge, is_twisted::Bool)
    e.is_twisted = is_twisted
end

trifaces(e::DualEdge) = (e.triangles[1],e.triangles[2])


struct TriFace
    id :: Int  # the unique index number of this face in the DeltaComplex
    points :: MVector{3, Int} #corners x,y,z
    edge_is_anticlockwise :: MVector{3, Bool}  #edge orientations xy, yz, zx  <- true if that order i.e. yx would be false
    edges :: Vector{DualEdge}  #edge xy, yz, zx

    TriFace(id::Int, x::Int, xy::Bool, y::Int, yz::Bool, z::Int, zx::Bool) = new(id, [x,y,z], [xy,yz,zx], Vector{DualEdge}(undef,3))
    TriFace(id::Int, x::Int, y::Int, z::Int) = new(id, [x,y,z], [true,true,true], Vector{DualEdge}(undef,3) )#MVector{3, DualEdge}(undef))
    TriFace(id::Int, points :: Vector{<:Integer}, edge_is_anticlockwise :: Vector{Bool}, edges:: Vector{DualEdge}) = new(id,points, edge_is_anticlockwise, edges)
end

has_point(T::TriFace, x::Integer) = (x in T.points) :: Bool
get_points(T::TriFace) = return Tuple(T.points) :: Tuple{Int, Int, Int}
get_point(T::TriFace, point_index::Integer) = T.points[point_index]
edge_is_anticlockwise(T::TriFace, edge_index::Integer) = T.edge_is_anticlockwise[edge_index] :: Bool
set_edge_anticlockwise!(T::TriFace, edge_index::Integer, is_anticlockwise::Bool) = T.edge_is_anticlockwise[edge_index] = is_anticlockwise :: Bool 
set_edge!(T::TriFace, edge_index::Integer, edge::DualEdge) = (T.edges[edge_index] = edge)
get_edge(T::TriFace, edge_index::Integer) = T.edges[edge_index] :: DualEdge
get_edges(T::TriFace) = T.edges :: Vector{DualEdge}
get_neighbor(T::TriFace, edge_index::Integer) = get_other_endpoint(T.edges[edge_index], T.id) :: Int

function Base.show(io::IO, mime::MIME"text/plain", T::TriFace)
    print(io, string("TriFace #",T.id, ": Points(", T.points[1]," ", T.points[2]," ", T.points[3],")   SideOrientations("))
    for b in T.edge_is_anticlockwise
        if b
            print(io, " ↺ ")
        else
            print(io, " ↻ ")
        end
    end
    print(io, string(")   Neighbours(",get_neighbor(T,1), " ",get_neighbor(T,2), " ",get_neighbor(T,3), ")"))
end

function get_triangle_edge(T::TriFace, edge_index::Integer)
    edge = (T.points[edge_index], T.points[(edge_index%3) + 1])
    if edge_is_anticlockwise(T, edge_index)
        return edge
    else
        return reverse(edge)
    end
end

"""
    function get_edge_index(T :: TriFace, e:: DualEdge)
returns the index 1,2 or 3 of the edge e in the TriFace T. If e is not incident to T this function returns -1.


function get_edge_index(T :: TriFace, e:: DualEdge)
    if T.v[1] in e.edge
        if T.v[2] in e.edge
            return 1
        elseif T.v[3] in e.edge
            return 3
        end
    elseif T.v[2] in e.edge && T.v[3] in e.edge
        return 2
    end
    return -1
end

function get_edge(T:: TriFace, i::Int) ::DualEdge
    return T.E[i]
end
"""

"""
    struct DeltaComplex
    
A Graph datastructure representing a triangulation of a surface. \\
Vertices are triangular faces. Every vertex has three edges incident to it.
"""
struct DeltaComplex
    V :: Array{TriFace, 1}
    E :: Array{DualEdge, 1}    
    num_points :: Base.RefValue{Int}
    #E :: Array{Tuple{Int, Int, Bool}, 1} #t1, t2, is_twisted_edge
    #is_twisted_edge :: Dict{Tuple{Int, Int}, Bool} #true if twisted, else false

    function DeltaComplex()
        new(TriFace[], DualEdge[], Ref(1))
    end
end

add_vertex!(D::DeltaComplex, v::TriFace) = push!(D.V, v)
get_vertex(D::DeltaComplex, i::Integer)::TriFace = D.V[i] 
get_vertices(D::DeltaComplex) = D.V ::Array{TriFace, 1}
add_edge!(D::DeltaComplex, e::DualEdge) = push!(D.E, e)
get_edge(D::DeltaComplex, i::Integer) = D.E[i]
get_edges(D::DeltaComplex) = D.E
remove_edge!(D::DeltaComplex, e::DualEdge) = remove!(D.E, e)
set_num_points!(D::DeltaComplex, num_points) = setindex!(D.num_points, num_points)

"""
    get_num_points(D::DeltaComplex)

Return the number of points in the triangulation defined by `D`
"""
get_num_points(D::DeltaComplex) = np(D) :: Int
np(D::DeltaComplex) = getindex(D.num_points) :: Int

"""
    get_num_trifaces(D::DeltaComplex)

Return the number of vertices(i.e. triangular faces) in `D`.
"""
get_num_trifaces(D::DeltaComplex) = nv(D)
nv(D::DeltaComplex) = length(D.V) :: Integer

"""
    get_num_edges(D::DeltaComplex)
"""
get_num_edges(D::DeltaComplex) = ne(D)
ne(D::DeltaComplex) = length(D.E) :: Integer

function Base.show(io::IO, mime::MIME"text/plain", D::DeltaComplex)
    if is_orientable(D)
        println(io, string("DeltaComplex on orientable surface of genus ", genus(D), " with ", np(D), " point",(if np==1 "s";else "";end)))
    else
        println(io, string("DeltaComplex on non-orientable surface of demigenus ", demigenus(D), " with ", np(D), " point", (if np==1 "s";else "";end)))
    end
    println(io, string(nv(D), " TriFaces:"))
    for T in get_vertices(D)
        print(io," "); show(io, mime, T); println(io)
    end
    println(io, string(ne(D), " DualEdges:"))
    for e in D.E
        print(io," "); show(io, mime, e); println(io)
    end  
end

"""
    euler_characteristic(D::DeltaComplex)

Compute the euler characteristic of the DeltaComplex `D`:\\
`` X = #vertices - #edges + #faces ``
"""
euler_characteristic(D::DeltaComplex) = np(D) - ne(D) + nv(D) :: Int

"""
    genus(D::DeltaComplex)

Compute the genus of the DeltaComplex `D` if D forms an orientable surface.
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

The demigenus (``k``) of a connected non-orientable closed surface is defined via the euler characteristic ()``X``) :\\
``k = 2 - X``
"""
function demigenus(D::DeltaComplex)
    if is_orientable(D)
        throw(ArgumentError("Cannot compute the demigenus of a orientable surface"))
    end
    return (2 - euler_characteristic(D))
end

"""
    adjacency_matrix(D::DeltaComplex) :: Matrix{<:Integer}

Compute the adjacency matrix of the delta complex `D`.
"""
function adjacency_matrix_deltaComplex(D::DeltaComplex) :: Matrix{<:Integer}
    A = zeros(Int, get_num_trifaces(D), get_num_trifaces(D))
    foreach(e-> (A[e.triangles[1], e.triangles[2]] = 1; A[e.triangles[2], e.triangles[1]] = 1), get_edges(D))
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
        ), get_vertices(D))
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
glue_faces_along_edge!(D::DeltaComplex, T1:: TriFace, edge_index_1::Integer, T2:: TriFace, edge_index_2::Integer) = glue_faces_along_edge!(D, T1.id, edge_index_1, T2.id, edge_index_2)
function glue_faces_along_edge!(D::DeltaComplex, t1:: Integer, edge_index_1::Integer, t2:: Integer, edge_index_2::Integer)
    T1 = get_vertex(D, t1)
    T2 = get_vertex(D, t2)

    #set_neighbor!(T1, edge_index_1, t2)
    #set_neighbor!(T2, edge_index_2, t1)

    x1,y1 = get_triangle_edge(T1, edge_index_1)
    x2,y2 = get_triangle_edge(T2, edge_index_2)
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

    e =  DualEdge(t1, edge_index_1, t2, edge_index_2, (edge_is_anticlockwise(T1, edge_index_1) == edge_is_anticlockwise(T2, edge_index_2)) )
    add_edge!(D,e)
    set_edge!(T1, edge_index_1, e)
    set_edge!(T2, edge_index_2, e)
      
end


function merge_points!(D::DeltaComplex, x::Integer, y::Integer)
    if x > y
        x, y = y, x
    end
    if x!=y
        set_num_points!(D, get_num_points(D)-1)
    end
    return rename_points!(D, y, x)
end

function rename_points!(D::DeltaComplex, x_old::Integer, x_new::Integer)
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
    for i = 2:num_points
        subdivide(D, i-1)
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
        elseif i<n
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
                    Delta-Complexes with boundary haven't been implemented yet.\n 
                Make sure, that every scaffolding edge is assigned to exactly 1 other edge."))
            else
                throw(ArgumentError(string(i)*" appears more than twice in s. \n
                    Each edge can only be assigned to 1 other edge"))
            end
        end
    end

    D = createDeltaScaffold(size(s,1))

    for a in s_abs_unique
        i,j = findall(σ -> σ==a , s_abs)
        Ti, edge_idx_i = get_TriFace_and_rel_edge(D, i)
        Tj, edge_idx_j = get_TriFace_and_rel_edge(D, j)

        if s[i] < 0
            set_edge_anticlockwise!(Ti, edge_idx_i, false)
        end
        if s[j] < 0
            set_edge_anticlockwise!(Tj, edge_idx_j, false)
        end
        
        glue_faces_along_edge!(D, Ti.id, edge_idx_i, Tj.id, edge_idx_j)
    end
    return D
end


function createDeltaScaffold(num_vertices :: Integer)
    if num_vertices%2 != 0 
        throw(ArgumentError("num_vertices has to be a multiple of 2"))
    elseif  num_vertices<=0
        throw(ArgumentError("num_vertices has to be strictly positive"))
    end

    n = num_vertices
    D = DeltaComplex()
    set_num_points!(D, n)

    T1 = TriFace(1, 1, true, 2, true, n, true)
    T2 = TriFace(2, 2, true, n-1, true, n, false)
    add_vertex!(D, T1)
    add_vertex!(D, T2)
    glue_faces_along_edge!(D, T1, 2, T2, 3)

    for i = 1 : (num_vertices-4) ÷ 2
        T1 = TriFace(2*i+1, i+1, true, i+2, true, n-i, false)
        add_vertex!(D, T1)
        glue_faces_along_edge!(D, T2, 1, T1, 3)
        T2 = TriFace(2*i+2, i+2, true, n-i-1, true, n-i, false)
        add_vertex!(D, T2)
        glue_faces_along_edge!(D, T1, 2, T2, 3)
    end

    return D
end


"""
    subdivide(D::DeltaComplex, t::Integer)

Add a point to the inside of the `t`-th TriFace and connect it to each corner.
"""
function subdivide(D::DeltaComplex, t::Integer)
    set_num_points!(D, get_num_points(D)+1)
    n = get_num_trifaces(D)
    T = get_vertex(D, t)    
    v1,v2,v3 = get_points(T)
    v0 = get_num_points(D)
    
    edges_out = filter(e -> (t in e.triangles ),  D.E)

    filter!(e -> !(t in e.triangles ),  D.E)
    T1 = TriFace(T.id, [v1, v2, v0], [edge_is_anticlockwise(T,1), true, false], Vector{DualEdge}(undef,3)) # [get_neighbor(T,1) , n+1, n+2])
    T2 = TriFace(n+1, [v2, v3, v0], [edge_is_anticlockwise(T,2), true, false], Vector{DualEdge}(undef,3)) #[get_neighbor(T,2) , n+2, T.id])
    T3 = TriFace(n+2, [v3, v1, v0], [edge_is_anticlockwise(T,3), true, false], Vector{DualEdge}(undef,3)) #[get_neighbor(T,3) , T.id, n+1])

    D.V[t] = T1
    push!(D.V, T2, T3)
  
    for i in eachindex(edges_out)
        if edges_out[i].triangles[1] != T.id
            edges_out[i] = reverse(edges_out[i])
        end
    end        
    sort!(edges_out, by = (e -> e.resp_sides[1]))

    if length(edges_out) == 2  # the triangle shares an edge with itself
        i = 1
        if edges_out[2].triangles[1] == edges_out[2].triangles[2]
            i = 2
        end
        push!(edges_out, reverse!(deepcopy(edges_out[i])))
    end
    #glue outer edges together
    glue_faces_along_edge!(D, T1.id, 1, get_neighbor(T,1), edges_out[1].resp_sides[2])
    glue_faces_along_edge!(D, T2.id, 1, get_neighbor(T,2), edges_out[2].resp_sides[2])
    glue_faces_along_edge!(D, T3.id, 1, get_neighbor(T,3), edges_out[3].resp_sides[2])

    #glue inner edges together
    glue_faces_along_edge!(D, T1.id, 2, T2.id, 3)
    glue_faces_along_edge!(D, T2.id, 2, T3.id, 3)
    glue_faces_along_edge!(D, T3.id, 2, T1.id, 3)
    return D
end

"""
    is_flippable(D::DeltaComplex, e_id::Integer)
    is_flippable(e::DualEdge)

Return true if the given edge is can be flipped.\\ 
This is only the case if the edge does not connect a triangle face to itself.
"""
is_flippable(D::DeltaComplex, e_id::Integer) = is_flippable(get_edge(D,e_id))
function is_flippable(e::DualEdge)
    return e.triangles[1] != e.triangles[2]
end

"""
    flip!(D::DeltaComplex, e_id::Integer)
    flip!(D::DeltaComplex, e::DualEdge)

Flip the given edge in D
"""
flip!(D::DeltaComplex, e_id::Integer) = flip!(D, get_edge(D,e_id))
function flip!(D::DeltaComplex, e::DualEdge)
    T1 = get_vertex(D, e.triangles[1])
    T2 = get_vertex(D, e.triangles[2])
    
    #find edges a,b,c,d and points x,y,z,q
    if e.resp_sides[1]==1
        t,a,b = get_edges(T1)
        z,x,y = get_points(T1)
        rot1 = 0
    elseif e.resp_sides[1]==2
        b,t,a = get_edges(T1)
        y,z,x = get_points(T1)
        rot1 = 1
    elseif e.resp_sides[1]==3
        a,b,t = get_edges(T1)
        x,y,z = get_points(T1)
        rot1 = -1
    end

    if e.resp_sides[2]==1
        t,c,d = get_edges(T2)
        q = get_point(T2,3)
        rot2 = 0
    elseif e.resp_sides[2]==2
        d,t,c = get_edges(T2)
        q = get_point(T2,1)
        rot2 = 1
    elseif e.resp_sides[2]==3
        c,d,t = get_edges(T2)
        q = get_point(T2,2)
        rot2 = -1
    end

    #The flip itself
    T1.edges .= e,b,c
    T1.points .= q,y,z
    T2.edges .= e,d,a
    T2.points .= y,q,x

    e.resp_sides.= 1,1
    update_endpoint!(a, T1.id, 2 + rot1, T2.id, 3)
    update_endpoint!(b, T1.id, (3+rot1-1)%3 + 1, T1.id, 2)
    update_endpoint!(c, T2.id, 2 + rot2, T1.id, 3)
    update_endpoint!(d, T2.id, (3+rot2-1)%3 + 1, T2.id, 2)
    #TODO either update side orientations or remove them alltogether

    return D
end

"""
    is_orientable(D::DeltaComplex)
"""
function is_orientable(D::DeltaComplex)
    color = zeros(Int8, nv(D))
    v0 = 1
    color[v0] = 1
    queue = [v0]
    while !isempty(queue)
        v = pop!(queue)
        for i = 1:3
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


#TODO check if everything works with twisted edges
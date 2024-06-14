"""
    struct FGVertexCandidate

A candidate for a possible vertex in a flip graph. 

It hasn't yet been decided if this candidate is a new vertex in the graph, or isomorph to an already existing `DeltaComplex`.
"""
struct FGVertexCandidate
    D::DeltaComplex

    num_point_perms :: Int
    num_edge_perms :: Int

    p_degrees :: Vector{Int}

    multi_adjacency_matrix_triangulation
    adjacency_matrix_deltacomplex

    function FGVertexCandidate(D::DeltaComplex, labeled_points::Bool)
        if !labeled_points
            point_perms = mcKay_points(D)
            rename_points!(D, point_perms[1])
            rename_vertices!(D, mcKay_vertices(D, only_one=true)[1])
            edge_perms = mcKay_edges(D)
            rename_edges!(D, edge_perms[1])
            new(D, length(point_perms), length(edge_perms), point_degrees(D),
            multi_adjacency_matrix_triangulation(D), adjacency_matrix_deltacomplex(D))
        else
            rename_vertices!(D, mcKay_vertices(D, only_one=true)[1])
            edge_perms = mcKay_edges(D)
            rename_edges!(D, edge_perms[1])
            new(D, 1, length(edge_perms), point_degrees(D),
            multi_adjacency_matrix_triangulation(D), adjacency_matrix_deltacomplex(D))
        end 
    end

    function FGVertexCandidate(D::DeltaComplex, labeled_points::Bool, edge_ids::Vector{<:Integer})
        if !labeled_points
            point_perms = mcKay_points(D)
            rename_points!(D, point_perms[1])
        else
            point_perms = [collect(1:np(D))]
        end
        rename_vertices!(D, mcKay_vertices(D,only_one=true)[1])
        edge_perms = mcKay_edges(D)
        rename_edges!(D, edge_perms[1])
        return new(D, length(point_perms), length(edge_perms), point_degrees(D),
            multi_adjacency_matrix_triangulation(D), adjacency_matrix_deltacomplex(D)),
            edge_perms[1][edge_ids]
    end
end

"""
    struct FGVertex

A *vertex* in a `FlipGraph`. 

A `FGVertex` is composed of a representant (`DeltaComplex`) of the isotopy class of that vertex. 
The representant has been relabeled with one of the canonical labeling obtained by an adaptation of *McKay's Algorithm*.
Additionally, the `FGVertex` contains the number of labeling that are output by the respective *McKay's Algorithms*.
"""
struct FGVertex
    D :: DeltaComplex
    id :: Int

    np :: Int
    nv :: Int
    ne :: Int

    point_perms :: Vector{Vector{Int}}
    vertex_perms :: Vector{Vector{Vector{Int}}}
    edge_perms :: Vector{Vector{Vector{Vector{Int}}}}

    p_degrees :: Vector{Int} #sorted if labeled_points=false

    multi_adjacency_matrix_triangulation :: Matrix{Int32}
    adjacency_matrix_deltacomplex :: Matrix{Int32}

    function FGVertex(D::DeltaComplex, id::Integer, labeled_points::Bool)
        A_deltacomplex = adjacency_matrix_deltacomplex(D)
        if labeled_points
            point_perms = [collect(1:np(D))]
        else
            point_perms = mcKay_points(D)
        end
        vertex_perms = Vector{Vector{Vector{Int}}}(undef, length(point_perms)) 
        edge_perms = [Vector{Vector{Vector{Int}}}() for i in 1:length(point_perms)]
        for i_p in eachindex(point_perms)
            if labeled_points
                vertex_perms[i_p] = mcKay_vertices(D, A_deltacomplex = A_deltacomplex)
            else
                vertex_perms[i_p] = mcKay_vertices(D, A_deltacomplex, point_perms[i_p])
            end
            for i_v in eachindex(vertex_perms[i_p])
                push!(edge_perms[i_p], mcKay_edges(D, point_perms[i_p], vertex_perms[i_p][i_v]))
            end
        end 
        new(D, id, np(D), nv(D), ne(D), point_perms, vertex_perms, edge_perms, 
        (labeled_points ? point_degrees(D) : sort!(point_degrees(D))), 
        multi_adjacency_matrix_triangulation(D), A_deltacomplex)
    end

    function FGVertex(fgvc::FGVertexCandidate, id::Integer, labeled_points::Bool)
        D = deepcopy(fgvc.D)
        if labeled_points
            point_perms = [collect(1:length(fgvc.p_degrees))]
        else
            point_perms = mcKay_points(D)
        end
        vertex_perms = Vector{Vector{Vector{Int}}}(undef, length(point_perms)) 
        edge_perms = [Vector{Vector{Vector{Int}}}() for i in eachindex(point_perms)]
        for i_p in eachindex(point_perms)
            if labeled_points
                vertex_perms[i_p] = mcKay_vertices(D, A_deltacomplex = fgvc.adjacency_matrix_deltacomplex)
            else
                vertex_perms[i_p] = mcKay_vertices(D, fgvc.adjacency_matrix_deltacomplex, point_perms[i_p])
            end
            for i_v in eachindex(vertex_perms[i_p])
                push!(edge_perms[i_p], mcKay_edges(D, point_perms[i_p], vertex_perms[i_p][i_v]))
            end
        end 
        new(D, id, np(D), nv(D), ne(D), point_perms, vertex_perms, edge_perms, 
        fgvc.p_degrees, fgvc.multi_adjacency_matrix_triangulation, fgvc.adjacency_matrix_deltacomplex)
    end
end

function get_point_perm(fgv::FGVertex, i_p::Integer)
    return fgv.point_perms[i_p]
end

function get_point_perms(fgv::FGVertex)
    return fgv.point_perms
end

function get_vertex_perm(fgv::FGVertex, i_p::Integer, i_v::Integer)
    return fgv.vertex_perms[i_p][i_v]
end

function get_vertex_perms(fgv::FGVertex, i_p::Integer)
    return fgv.vertex_perms[i_p]
end

function get_edge_perms(fgv::FGVertex, i_p::Integer, i_v::Integer)
    return fgv.edge_perms[i_p][i_v]
end

function get_edge_perm(fgv::FGVertex, i_p::Integer, i_v::Integer, i_e::Integer)
    return fgv.edge_perms[i_p][i_v][i_e]
end

#function get_D(fgv::FGVertex, i_p::Integer)
#    D = deepcopy(fgv.D)
#    return rename_points!(D, fgv.point_perms[i_p])
#end
#
#function get_D(fgv::FGVertex, i_p::Integer, i_v::Integer)
#    D = deepcopy(fgv.D)
#    rename_points!(D, fgv.point_perms[i_p])
#    return rename_vertices!(D, fgv.vertex_perms[i_p][i_v])
#end
#
#function get_D(fgv::FGVertex, i_p::Integer, i_v::Integer, i_e::Integer)
#    D = deepcopy(fgv.D)
#    rename_points!(D, fgv.point_perms[i_p])
#    rename_vertices!(D, fgv.vertex_perms[i_p][i_v])
#    return rename_edges!(D, fgv.edge_perms[i_p][i_v][i_e])
#end


"""
    struct FlipGraph <: AbstractGraph{Int}

A *graph* representing the **flip graph** of a **Δ-Complex**.

Vertices are isotopy classes of triangulations of the same surface.\\
Two vertices are linked by an edge, if the respective triangulations differ only by a single flip.
"""
struct FlipGraph <: AbstractGraph{Int32}    
    V::Vector{FGVertex}
    adjList::Vector{Vector{Int32}}
    labeled_points::Bool

    function FlipGraph(labeled_points::Bool)
        new(Vector{FGVertex}(), Vector{Vector{Int32}}(), labeled_points)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", G::FlipGraph)
    print(io, string("modular FlipGraph with ", nv(G) , " vertices and ", ne(G), " edges")); 
end

"""
    edges(G::FlipGraph) :: Vector{Edge}

Construct a list containing all the edges in `G`.
"""
function edges(G::FlipGraph) ::Vector{Edge}
    return collect(Edge(Int32(i),j) for i in eachindex(G.V) for j in G.adjList[i] if i<j)
end 

edgetype(G::FlipGraph) = SimpleEdge{Int32}

"""
    has_edge(G::FlipGraph, e::Edge) :: Bool

Return `true` if `e` is an edge in `G`.
"""
has_edge(G::FlipGraph, e::Edge)::Bool = (dst(e) ∈ G.adjList[src(e)])

"""
    has_edge(G::FlipGraph, s::Integer, d::Integer) :: Bool

Return `true` if there is an edge between the `s`-th and `d`-th vertex in `G`.
"""
has_edge(G::FlipGraph, s::Integer, d::Integer) :: Bool = (d ∈ G.adjList[s])

"""
    has_edge(G::FlipGraph, v1::FGVertex, v2::FGVertex) :: Bool

Return `true` if there is an edge between `v1` and `v2` in `G`.
"""
has_edge(G::FlipGraph, v1::FGVertex, v2::FGVertex) :: Bool = (v2.id ∈ G.adjList[v1.id])


"""
    has_vertex(G::FlipGraph, v::Integer) :: Bool

Return true if `v` is a valid index of a vertex in `G`.
"""
has_vertex(G::FlipGraph, v::Integer) :: Bool = (1 <= v <= nv(G))

"""
    neighbors(G::FlipGraph, v::Integer) :: Vector{Int32}

Return a list of the indices of all the neighboring vertices of the `v`-th vertex.
"""
neighbors(G::FlipGraph, v::Integer) :: Vector{Int32} = G.adjList[v]
inneighbors(G::FlipGraph, v) = neighbors(G,v)
outneighbors(G::FlipGraph, v) = neighbors(G,v)

"""
    ne(G::FlipGraph) :: Int

Return the number of edges in `G`.
"""
ne(G::FlipGraph) :: Int = sum(size(G.adjList[i],1) for i in eachindex(G.adjList))÷2

"""
    nv(G::FlipGraph) :: Int

Return the number of vertices in `G`.
"""
nv(G::FlipGraph) :: Int = length(G.V)

"""
    vertices(G::FlipGraph) :: Vector{FGVertex}

Return the list of all the vertices that have been constructed in `G`.
"""
vertices(G::FlipGraph) :: Vector{FGVertex} = G.V

"""
    get_vertex(G::FlipGraph, id::Integer) :: Vector{FGVertex}

Return the `id`-th vertex the flip graph `G`.
"""
get_vertex(G::FlipGraph, id::Integer) :: FGVertex = G.V[id]

"""
    vertices_deltacomplex(G::FlipGraph) :: Vector{DeltaComplex}

Construct a list of all the `DeltaComplex`s which form the vertices in `G`.
"""
vertices_deltacomplex(G::FlipGraph) :: Vector{DeltaComplex} = collect(DeltaComplex, v.D for v in vertices(G))

is_directed(G::FlipGraph) ::Bool = false
is_directed(::Type{FlipGraph}) ::Bool = false


function add_edge!(G::FlipGraph, v, w) 
    if !has_edge(G, v, w) && v!=w
        push!(G.adjList[v], w)
        push!(G.adjList[w], v)
    end
end

function add_vertex!(G::FlipGraph, fgv::FGVertex) :: FGVertex
    push!(G.V, fgv)
    push!(G.adjList,[])
    return fgv
end

function add_vertex!(G::FlipGraph, D::DeltaComplex) :: FGVertex
    fgv = FGVertex(D, nv(G)+1, G.labeled_points)
    return add_vertex!(G, fgv)
end

function add_vertex!(G::FlipGraph, fgvc::FGVertexCandidate) :: FGVertex
    fgv = FGVertex(fgvc, nv(G)+1, G.labeled_points)
    return add_vertex!(G, fgv)
end

#function build_vertex_tree(nps::Integer, nes::Integer, labeled_points::Bool)
#    if nps == 1
#        return Vector{FGVertex}()
#    end
#    T = Vector{Any}(collect(2*nes-1:-1:(nps-1)))
#    t = T
#    d = 1
#    while d < nps
#        for i in eachindex(t)
#            if d < nps-1
#                t[i] = Vector{Any}(collect(t[i]-1):-1:(nps-d))
#            else
#                t[i] = Vector{Any}(undef, t[i])
#            end
#        end
#        d += 1
#        t
#    end
#end


"""
    flipgraph_modular(g::Integer, p::Integer; kwargs..) :: FlipGraph
    
Construct the **modular flip graph** for a genus `g` closed orientable surface with `p` points on it.  

# Arguments
- `labeled_points :: Bool = true` : If set to `false`, the isomorphism also includes a renaming of the points. 
"""
function flipgraph_modular(g::Integer, p::Integer; labeled_points::Bool=true)
    return flipgraph_modular(deltacomplex(g, p), labeled_points=labeled_points)
end

"""
    flipgraph_modular(D::DeltaComplex; kwargs..)
    
Construct the **modular flip graph** for the closed orientable surface defined by the *Δ-complex* `D`.  

# Arguments
- `labeled_points :: Bool = true` : If is set to `false`, then the isomorphism also includes a renaming of the points. 
- `depth :: Integer = ∞` : Determines the depth to which the flip grap should be constructed. i.e. up to which distance from `D`. 
"""
function flipgraph_modular(D::DeltaComplex; depth::Integer = typemax(Int), labeled_points::Bool = true)
    G = FlipGraph(labeled_points)
    nes = ne(D); nps = np(D)#TODO add the first vertex to the tree    
    fgv_first = add_vertex!(G, D)

    #rooted tree whose i-th nodes branch-off according to the possible degrees of the i-th point. 
    #The leafs contain all the vertices that have the same point_degrees
    if nps > 1
        if labeled_points
            vertex_tree = Vector{Any}(undef, 2*nes-(nps-1)) 
        else
            vertex_tree = Vector{Any}(undef, (2*nes)÷nps) 
        end
    else
        vertex_tree = Vector{FGVertex}([fgv_first])
    end

    vtn = vertex_tree
    for i in 1:nps-1 #fgvc.p_degrees[1:end-1] 
        if i != nps-1
            if labeled_points
                vtn[fgv_first.p_degrees[i]] = Vector{Any}(undef, 2*nes - sum(fgv_first.p_degrees[1:i])-(nps-i-1))
            else
                vtn[fgv_first.p_degrees[i]] = Vector{Any}(undef, (2*nes - sum(fgv_first.p_degrees[1:i]))÷(nps-i))
            end
        else
            vtn[fgv_first.p_degrees[i]] = Vector{FGVertex}([fgv_first])
        end
        vtn = vtn[fgv_first.p_degrees[i]]
    end

    queue = Vector{Tuple{FGVertex, Int, Int}}()  #(D, fgv, index, depth)
    push!(queue, (fgv_first, 1, 0))
    while !isempty(queue)
        fgv_first, ind_D, d = popfirst!(queue)
        D = deepcopy(fgv_first.D)
        edge_ids = collect(1:ne(D))
        for e in 1:ne(D)
            if is_flippable(D, edge_ids[e])
                flip!(D, edge_ids[e])
                fgvc, edge_ids = FGVertexCandidate(D, labeled_points, edge_ids)
                is_new = false
                vtn = vertex_tree
                for i in 1:nps-1 #fgvc.p_degrees[1:end-1] 
                    if !isassigned(vtn, fgvc.p_degrees[i])
                        if i != nps-1
                            if labeled_points
                                vtn[fgvc.p_degrees[i]] = Vector{Any}(undef, 2*nes - sum(fgvc.p_degrees[1:i])-(nps-i-1))
                            else
                                vtn[fgvc.p_degrees[i]] = Vector{Any}(undef, (2*nes - sum(fgvc.p_degrees[1:i]))÷(nps-i))
                            end
                        else
                            vtn[fgvc.p_degrees[i]] = Vector{FGVertex}()
                        end
                        is_new = true
                    end
                    vtn = vtn[fgvc.p_degrees[i]]
                end
                if !is_new
                    is_new = true
                    for v in vtn
                        if is_isomorphic(fgvc, v; labeled_points=labeled_points)
                            add_edge!(G, ind_D, v.id)
                            is_new = false
                            break
                        end
                    end
                end
                if is_new #add the new DeltaComplex to the Graph
                    add_vertex!(G, fgvc)
                    vtn = vertex_tree
                    for i in fgvc.p_degrees[1:end-1]
                        vtn = vtn[i]
                    end
                    push!(vtn, G.V[end])
                    add_edge!(G, ind_D, length(G.V))
                    if d + 1 < depth
                        push!(queue, (G.V[end], length(G.V), d+1))
                    end
                end
                flip!(D, edge_ids[e])
            end
        end
    end
    return G
end


"""
    is_isomorphic(D1::DeltaComplex, D2::DeltaComplex; kwargs..) :: Bool

Return `true` if `D1` is isomorph to `D2` up to a renaming of the vertices, edges and if `labeled_points=false` also points.
    
# Arguments
- `labeled_points :: Bool = true` : If labeled_points is set to `false`, then the isomorphism also includes a renaming of the points. 
"""
function is_isomorphic(D1::DeltaComplex, D2::DeltaComplex; labeled_points::Bool = true) :: Bool
    nv(D1) == nv(D2) && ne(D1) == ne(D2) && np(D1) == np(D2) || return false
    D = deepcopy(D1)
    fgvc = FGVertexCandidate(D, labeled_points)
    fgv = FGVertex(D2, 0, labeled_points)
    return is_isomorphic(fgvc, fgv, labeled_points=labeled_points)
end


"""
    is_isomorphic(candidate::FGVertexCandidate, fgv::FGVertex; kwargs...) -> Bool

Return `true` if `candidate` is in the isotopy class of `fgv`. 

# Arguments
- `labeled_points :: Bool = true` : If is set to `false`, then the isomorphism would also allow a relabeling of the points. 
"""
function is_isomorphic(candidate::FGVertexCandidate, fgv::FGVertex; labeled_points::Bool = true) :: Bool
    if candidate.p_degrees != fgv.p_degrees#TODO, do planar idea with only checking the right candidates
        return false
    end
    #point mapping
    if !labeled_points
        permutations_points = get_point_perms(fgv)
        valid_permutation_points = trues(length(permutations_points))
        if length(permutations_points) != candidate.num_point_perms
            return false
        end
        for i in eachindex(permutations_points)
            p = permutations_points[i]
            if !matrix_equal(fgv.multi_adjacency_matrix_triangulation, candidate.multi_adjacency_matrix_triangulation, p) 
                valid_permutation_points[i] = false
            end
        end
        if !any(valid_permutation_points) 
            return false
        end
    else
        if !matrix_equal(fgv.multi_adjacency_matrix_triangulation, candidate.multi_adjacency_matrix_triangulation)
            return false
        end
        permutations_points = [collect(1:fgv.np)]
        valid_permutation_points = [true]
    end

    for i_p in eachindex(permutations_points)
        #Triface mapping
        permutations_trifaces = get_vertex_perms(fgv, i_p)
        valid_permutation_vertices = trues(length(permutations_trifaces))
        for i in eachindex(permutations_trifaces) 
            p = permutations_trifaces[i]
            if !matrix_equal(fgv.adjacency_matrix_deltacomplex, candidate.adjacency_matrix_deltacomplex, p)
                valid_permutation_vertices[i] = false
            end
        end
        if !any(valid_permutation_vertices)
            continue
        end

        for i_v in eachindex(permutations_trifaces)
            if !valid_permutation_vertices[i_v]
                continue
            end

            #DualEdge mapping
            for i_e in eachindex(get_edge_perms(fgv,i_p,i_v))
                if is_all_equivalent(candidate.D, fgv, get_point_perm(fgv, i_p), get_vertex_perm(fgv, i_p, i_v), get_edge_perm(fgv, i_p, i_v, i_e))                 
                    return true
                end
            end
        end
    end
    return false
end
 
"""
    is_all_equivalent(D::DeltaComplex, fgv::FGVertex, p_p::Vector{T}, p_v::Vector{T}, p_e::Vector{T}) where T<:Integer

Return `true` if `D` is identical to the `DeltaComplex` defined by the vertex fgv, whose `DeltaComplex` has been renamed by the `p_p`-th point permutation, `p_v`-th vertex permutation, `p_e`-th edge permutation.
"""
function is_all_equivalent(D::DeltaComplex, fgv::FGVertex, p_p::Vector{T}, p_v::Vector{T}, p_e::Vector{T}) where T<:Integer
    for n in eachindex(D.V)
        T1 = D.V[p_v[n]]
        T2 = fgv.D.V[n]
        bo = true
        for offset in 0:2
            bo = true
            for i in 1:3
                j = (i+offset-1)%3 + 1
                if T1.points[i] != p_p[T2.points[j]]
                    bo = false
                    break
                elseif T1.edges[i].id != p_e[T2.edges[j].id]  #get_edge_id(T1, i) != p_e[get_edge_id(T2, j)]
                    bo = false
                    break
                end
            end
            if bo 
                break
            end
        end 
        if !bo
            return false
        end
    end
    return true
end


"""
    diameter(G::FlipGraph)

Compute the diameter of the `FlipGraph` `G`.
"""
function diameter(G::FlipGraph)
    return diameter(adjacency_matrix(G.adjList))
end

"""
    split(V::Vector{T}, degs::Vector{<:Integer}) :: Vector{Vector{T}} where T<:Integer

Split `V` into partitions according to their degrees from smallest to biggest.
"""
function split(V::Vector{T}, degs::Vector{<:Integer}) :: Vector{Vector{T}} where T<:Integer
    sV = Vector{Vector{T}}()
    unique_degs = sort!(unique(degs))
    j = 1 
    k = length(V)
    while k > 0
        W = Vector{T}()
        for i in eachindex(degs)
            if degs[i] == unique_degs[j]
                push!(W, V[i])
                k -= 1
            end
        end
        if !isempty(W)
            push!(sV, W)
        end
        j += 1
    end
    return sV
end


"""
    makeEquitable!(p::Vector{Vector{T}}, A::Matrix{<:Integer}) where T<:Integer
    
Replace partitions by sub-partitions as long as there are at least 2 elements in the same partition that may be differentiated by their relative degrees to another partition.
"""
function makeEquitable!(p::Vector{Vector{T}}, A::Matrix{<:Integer}) where T<:Integer
    i = 1; j = 1
    while i <= length(p)
        bo = true
        rd1 = relative_degree(A, p[i][1], p[j])
        for k in 2:length(p[i])
            if relative_degree(A, p[i][k], p[j]) != rd1
                bo = false
                break
            end
        end
        if !bo
            rDegs = relative_degrees(A, p[i], p[j])
            newVs = split(p[i], rDegs)        
            #replace the old partition by the new ones
            popat!(p,i)
            for V in reverse(newVs)
                insert!(p, i, V)
            end
            j = i
            i = 1
        else 
            j += 1
            if j > length(p)
                j = 1
                i += 1
            end
        end
    end 
    return p   
end

"""
    force_makeEquitable!(p::Vector{Vector{T}}, A::Matrix{<:Integer}) where T<:Integer

If it is no longer possible to makeEquitable, the first non-trivial set gets split by removing one of the elements and putting it in a singleton.
This is done with all the elements from said set, resulting in a bunch of possibilities. 
Each of the partitions is now again tried to be made equitable.
"""
function force_makeEquitable!(p::Vector{Vector{T}}, A::Matrix{<:Integer}, n::Integer) where T<:Integer
    queue = Vector{Vector{Vector{T}}}([p])
    leafs = Vector{Vector{T}}()
    while !isempty(queue)
        p = popfirst!(queue)
        i = 1
        while length(p[i]) == 1; i += 1 end #i = index of first partition that is not trivial
        for j in eachindex(p[i]) #replace the i-th partition by isolating its j-th element 
            pp = deepcopy(p)
            V = popat!(pp, i)
            insert!(pp, i, [V[j]])
            popat!(V, j)
            insert!(pp, i+1, V)
            makeEquitable!(pp, A)
            if length(pp) != n
                push!(queue, pp)
            else
                push!(leafs, invert_permutation(reduce(vcat, pp)))
            end
        end
    end
    return leafs
end

"""
    mcKay_points(D::DeltaComplex; kwargs..) :: Vector{Vector{Int32}}

Apply a version of McKay's canonical graph labeling algorithm to determine all possible permutations 
of the points which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that point 1 becomes point `p[1]`, point 2 becomes point `p[2]`,...

# Arguments
- `only_one :: Bool = false`, If set to true, the algorithm will stop after finding a single valid permutation.
"""
function mcKay_points(D::DeltaComplex; only_one::Bool=false) :: Vector{Vector{Int32}}
    A = multi_adjacency_matrix_triangulation(D)

    n = np(D)
    p = split(collect(Int32,1:n), degrees(A))
    makeEquitable!(p, A)

    if length(p) == n #there is only one canonical permutation
        return Vector{Vector{Int32}}([invert_permutation(reduce(vcat, p))])
    end

    while only_one     
        i = 1
        while length(p[i]) == 1 #i = index of first partition that is not trivial
            i += 1 
        end 
        V = popat!(p, i)
        insert!(p, i, [popat!(V, 1)])
        insert!(p, i+1, V)
        makeEquitable!(p, A)
        if length(p) == n
            return Vector{Vector{Int32}}([invert_permutation(reduce(vcat, p))])
        end
    end

    return force_makeEquitable!(p, A, n)  #Permutations
end

"""
    mcKay_vertices(D::DeltaComplex; , A_deltacomplex::Matrix{<:Integer}, point_perm::Vector{<:Integer}) :: Vector{Vector{Int32}}

Apply a version of McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the `TriFace`s which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that `TriFace` 1 becomes `TriFace` `p[1]`, `TriFace` 2 becomes `TriFace` `p[2]`, ...\\
If `only_one=true`, the algorithm stops after finding one valid permutation.

The vectors `point_perm` determins how the points of `D` are implied to have been renamed, without actually having been changed.
"""
function mcKay_vertices(D::DeltaComplex, A_deltacomplex::Matrix{<:Integer}, point_perm::Vector{<:Integer}) :: Vector{Vector{Int32}}
    np1 = np(D)
    np2 = np1*np1

    #give each combination of three points a distinct value (1,1,1)<(1,1,2)=(1,2,1)=(2,1,1)<(1,1,3)...  and (1,2,3)<(1,3,2)
    function pointValue(T::TriFace)
        p1,p2,p3 = point_perm[T.points].-1
        return min(p1*np2 + p2*np1 + p3, p2*np2 + p3*np1 + p1, p3*np2 + p1*np1 + p2)
    end

    n = nv(D)
    p = split(collect(Int32,1:n), collect(pointValue(T) for T in D.V))
    makeEquitable!(p,A_deltacomplex)
    if length(p) == n #there is only one canonical permutation
        return Vector{Vector{Int32}}([invert_permutation(reduce(vcat, p))])
    end
    #split the first partition that has more than 2 elements 
    return force_makeEquitable!(p, A_deltacomplex, n)
end


"""
    mcKay_vertices(D::DeltaComplex; kwargs..) :: Vector{Vector{Int32}}

Apply a version of McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the `TriFace`s which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that `TriFace` 1 becomes `TriFace` `p[1]`, `TriFace` 2 becomes `TriFace` `p[2]`, ...\\
If `only_one=true`, the algorithm stops after finding one valid permutation.

# Arguments
- `only_one :: Bool = false`, If set to true, the algorithm will stop after finding a single valid permutation.
- A\\_deltacomplex :: Matrix{<:Integer} = Matrix{Int32}(adjacency_matrix_deltacomplex(D)). If provided, the algorithm will use the adjacency matrix of the `DeltaComplex` `D`. If not, the algorithm would have to compute it itself.
"""
function mcKay_vertices(D::DeltaComplex; only_one::Bool=false, A_deltacomplex::Matrix{<:Integer} = Matrix{Int32}(adjacency_matrix_deltacomplex(D))) ::Vector{Vector{Int32}}
    np1 = np(D)
    np2 = np1*np1

    #give each combination of three points a distinct value (1,1,1)<(1,1,2)=(1,2,1)=(2,1,1)<(1,1,3)...  and (1,2,3)<(1,3,2)
    function pointValue(T::TriFace)
        p1,p2,p3 = T.points.-1
        return min(p1*np2 + p2*np1 + p3, p2*np2 + p3*np1 + p1, p3*np2 + p1*np1 + p2)
    end

    n = nv(D)
    p = split(collect(Int32, 1:n), collect(pointValue(T) for T in D.V))
    makeEquitable!(p, A_deltacomplex)

    if length(p) == n #|| only_one #there is only one canonical permutation
        return Vector{Vector{Int32}}([invert_permutation(reduce(vcat, p))])
    end
    
    while only_one     
        i = 1
        while length(p[i]) == 1 #i = index of first partition that is not trivial
            i += 1 
        end 
        V = popat!(p, i)
        insert!(p, i, [popat!(V, 1)])
        insert!(p, i+1, V)
        makeEquitable!(p, A_deltacomplex)
        if length(p) == n
            return Vector{Vector{Int32}}([invert_permutation(reduce(vcat, p))])
        end
    end
    
    return force_makeEquitable!(p, A_deltacomplex, n)
end

"""
    mcKay_edges(D::DeltaComplex, point_perm::Vector{<:Integer}, vertex_perm::Vector{<:Integer}) :: Vector{Vector{Int32}}

Apply McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the `DualEdge`s which give a canonical isomorphism class representant.

The vectors `point_perm` and `vertex_perm` respectively determin how the points and vertices of `D` are implied to have been renamed, without actually having been changed.

Return a vector of permutation vectors `p` such that `DualEdge` 1 becomes `DualEdge` `p[1]`, `DualEdge` 2 becomes `DualEdge` `p[2]`, ...\\
"""
function mcKay_edges(D::DeltaComplex, point_perm::Vector{<:Integer}, vertex_perm::Vector{<:Integer})
    b1 = max(np(D), nv(D))
    b2 = b1*b1

    #give each combination of two vertices and 2 points a distinct value 
    function edgeValue(d::DualEdge)
        t1,t2 = vertex_perm[d.triangles]
        p1,p2 = points(D, d)
        p1,p2 = point_perm[p1], point_perm[p2]
        return min(t1+t2*b1, t2+t1*b1)*b2 + min(p1+p2*b1, p2+p1*b1)
    end

    m = ne(D)
    p = split(collect(Int32, 1:m), collect(edgeValue(d) for d in D.E))
    i = 1
    while i <= length(p)
        if length(p[i]) > 1
            pp = popat!(p,i)
            insert!(p, i, pp)
            i += 1
        else
            i += 1
        end
    end
    #edges may only be in the same partition, if they share the same 2 endpoints and triangles
    #this is only possible for two edges at a time, they cannot be distinguished

    if length(p) == m #there is only one canonical permutation
        return Vector{Vector{Int32}}([invert_permutation(reduce(vcat, p))])
    end
    
    #split the first partition that has more than 2 elements 
    queue = Vector{Vector{Vector{Int32}}}([p])
    leafs = Vector{Vector{Int32}}()
    while !isempty(queue)
        p = popfirst!(queue)
        i = 1
        while length(p[i]) == 1; i += 1 end #i = index of first partition that is not trivial
        for j in eachindex(p[i]) #replace the i-th partition by isolating its j-th element 
            pp = deepcopy(p)
            V = popat!(pp, i)
            insert!(pp, i, [popat!(V, j)])
            insert!(pp, i+1, V)
            if length(pp) != m
                push!(queue, pp)
            else
                push!(leafs, invert_permutation(reduce(vcat, pp)))
            end
        end
    end
    return leafs   #sigma_pi' 
end

"""
    mcKay_edges(D::DeltaComplex; kwargs..) :: Vector{Vector{Int32}}

Apply McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the `DualEdge`s which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that `DualEdge` 1 becomes `DualEdge` `p[1]`, `DualEdge` 2 becomes `DualEdge` `p[2]`, ...\\

# Arguments
- `only_one :: Bool = false`, If set to true, the algorithm will stop after finding a single valid permutation.
"""
function mcKay_edges(D::DeltaComplex; only_one::Bool=false) :: Vector{Vector{Int32}}
    b1 = max(np(D), nv(D))
    b2 = b1*b1

    #give each combination of two vertices and 2 points a distinct value 
    function edgeValue(d::DualEdge)
        t1,t2 = d.triangles
        p1,p2 = points(D, d)
        return min(t1+t2*b1, t2+t1*b1)*b2 + min(p1+p2*b1, p2+p1*b1)
    end

    m = ne(D)
    p = split(collect(Int32, 1:m), collect(edgeValue(d) for d in D.E))
    i = 1
    while i <= length(p)
        if length(p[i]) > 1
            pp = popat!(p,i)
            insert!(p, i, pp)
            i += 1
        else
            i += 1
        end
    end
    #edges may only be in the same partition, if they share the same 2 endpoints and triangles
    #this is only possible for two edges at a time, they cannot be distinguished

    if length(p) == m #there is only one canonical permutation
        return Vector{Vector{Int32}}([invert_permutation(reduce(vcat, p))])
    end
    
    #split the first partition that has more than 2 elements 
    queue = Vector{Vector{Vector{Int32}}}([p])
    leafs = Vector{Vector{Int32}}()
    while !isempty(queue)
        p = popfirst!(queue)
        i = 1
        while length(p[i]) == 1; i += 1 end #i = index of first partition that is not trivial
        for j in eachindex(p[i]) #replace the i-th partition by isolating its j-th element 
            pp = deepcopy(p)
            V = popat!(pp, i)
            insert!(pp, i, [popat!(V, j)])
            insert!(pp, i+1, V)
            if length(pp) != m
                push!(queue, pp)
            else
                if only_one
                    return [invert_permutation(reduce(vcat, pp))]
                end
                push!(leafs, invert_permutation(reduce(vcat, pp)))
            end
        end
    end
    return leafs   #sigma_pi's
end
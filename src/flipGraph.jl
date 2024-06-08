"""
    struct FGVertex

A vertex in a flip graph. 

An `FGVertex` is composed of a representant(`DeltaComplex`) of the isotopy class of that vertex. 
The representant has been relabeld with one of the canonical labelings obtained by McKay's Algorithm.
Addidtionally, the `FGVertex` contains the number of labelings that are output by the respective McKay's Algorithms.
"""
struct FGVertex
    D::DeltaComplex

    num_point_perms :: Int
    num_vertex_perms :: Int
    num_edge_perms :: Int

    function FGVertex(D::DeltaComplex, fix_points::Bool)
        point_perms = mcKay_points(D)
        if !fix_points
            rename_points!(D, point_perms[1])
        end
        vertex_perms = mcKay_vertices(D)
        rename_vertices!(D, vertex_perms[1])
        edge_perms = mcKay_edges(D)
        rename_edges!(D, edge_perms[1])
        new(D, length(point_perms), length(vertex_perms), length(edge_perms))
    end
end

"""
    struct FlipGraph <: AbstractGraph{Int}

A Graph representing the flipgraph of a Δ-Complex.

Vertices are different triangulations of the same surface.\\
Two vertices are linked by an edge, if the respective graphs differ only by a single flip.
"""
struct FlipGraph <: AbstractGraph{Int}    
    V::Vector{FGVertex}
    adjList::Vector{Vector{Int}}
    fix_points::Bool

    function FlipGraph(fix_points::Bool)
        new(Vector{FGVertex}(), Vector{Vector{Int}}(), fix_points)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", G::FlipGraph)
    print(io, string("FlipGraph with ", nv(G) , " vertices and ", ne(G), " edges")); 
end

"""
    edges(G::FlipGraph) ::Vector{Edge}

Construct an array containing all the edges in `G`.
"""
function edges(G::FlipGraph) ::Vector{Edge}
    E = collect(Edge(i,j) for i in eachindex(G.V) for j in G.adjList[i])
    return filter!(e -> (src(e) > dst(e)), E)
end 

edgetype(G::FlipGraph) = SimpleEdge{Int}

"""
    has_edge(G::FlipGraph, e::Edge) -> Bool
"""
has_edge(G::FlipGraph, e::Edge)::Bool = (dst(e) ∈ G.adjList[src(e)])

"""
    has_edge(G::FlipGraph, s::Integer, d::Integer) -> Bool
"""
has_edge(G::FlipGraph, s::Integer, d::Integer)::Bool = (d ∈ G.adjList[s])

"""
    has_edge(G::FlipGraph, D1::DeltaComplex, D2::DeltaComplex) -> Bool
"""
function has_edge(G::FlipGraph, D1::DeltaComplex, D2::DeltaComplex)::Bool 
    return has_edge(G, findfirst(x->x==D1, G.V), findfirst(x->x==D2, G.V))
end

"""
    has_vertex(G::FlipGraph, v::Integer) -> Bool

Return true if `v` is a valid index of a vertex in `G`.
"""
has_vertex(G::FlipGraph, v::Integer) = (1 <= v <= nv(G))
inneighbors(G::FlipGraph, v) = G.adjList[v]
outneighbors(G::FlipGraph, v) = G.adjList[v]

"""
    neighbors(G::FlipGraph, v::Integer) -> Vector{Int}

Return a list of the indices of all the neighboring vertices of the `v`-th vertex.
"""
neighbors(G::FlipGraph, v::Integer) ::Vector{Int} = G.adjList[v]

"""
    ne(G::FlipGraph) -> Int

Return the number of edges in `G`.
"""
ne(G::FlipGraph) ::Int = sum(size(G.adjList[i],1) for i in eachindex(G.adjList))÷2

"""
    nv(G::FlipGraph) -> Int

Return the number of vertices in `G`.
"""
nv(G::FlipGraph) ::Int = length(G.V)

"""
    vertices(G::FlipGraph) -> Vector{DeltaComplex}

Return a list of all the vertices that have been constructed in `G`.
"""
vertices(G::FlipGraph) :: Vector{DeltaComplex} = G.V
is_directed(G::FlipGraph) = false
is_directed(::Type{FlipGraph}) = false


function add_edge!(G::FlipGraph, v, w) 
    if !has_edge(G, v, w) && v!=w
        push!(G.adjList[v],w)
        push!(G.adjList[w],v)
    end
end

function add_vertex!(G::FlipGraph, D::DeltaComplex) 
    fgv = FGVertex(D, G.fix_points)
    push!(G.V, fgv)
    push!(G.adjList,[])
    return fgv
end

function remove_edge!(G::FlipGraph, e::Edge)
    deleteat!(G.adjList[src(e)], findfirst(x->x==dst(e), G.adjList[src(e)]))
    deleteat!(G.adjList[dst(e)], findfirst(x->x==src(e), G.adjList[dst(e)]))
end


"""
    flipgraph_modular(D::DeltaComplex, depth::Integer; kwargs)
    
Construct the **FlipGraph** for the DeltaComplex `D`.  

# Arguments
- `fix_points::Bool=true` : If is set to `false`, then the isomorphism also includes a renaming of the points. 
"""
function flipgraph_modular(D::DeltaComplex, depth::Integer; fix_points::Bool = true)
    G = FlipGraph(fix_points)
    if !fix_points
        D = rename_points!(D, mcKay_points(D, only_one=true)[1])
    end
    fgv = add_vertex!(G, D)
    queue = Vector{Tuple{FGVertex, Int, Int}}()  #(D, index, depth)
    push!(queue, (fgv, 1, 0))
    while !isempty(queue)
        fgv, ind_D, d = popfirst!(queue)
        D = fgv.D
        for e in 1:ne(D)
            if is_flippable(D, e)
                new_D = flip(D, e)
                is_new = true
                for i in eachindex(G.V)
                    if is_isomorphic(G.V[i], new_D; fix_points=fix_points)
                        add_edge!(G, ind_D, i)
                        is_new = false
                        break
                    end
                end
                if is_new #add the new DeltaComplex to the Graph
                    add_vertex!(G, new_D)
                    add_edge!(G, ind_D, nv(G))
                    if d + 1 < depth
                        push!(queue, (G.V[end], nv(G), d+1))
                    end
                end
            end
        end
    end
    return G
end


"""
    is_isomorphic(D1::DeltaComplex, D2::DeltaComplex; kwargs) -> Bool

Return `true` if `D1` is isomorph to `D2` up to a renaming of the vertices, edges and if `fix_points=false` also points.
    
# Arguments
- `fix_points::Bool=true` : If fix_points is set to `false`, then the isomorphism also includes a renaming of the points. 
"""
function is_isomorphic(D1::DeltaComplex, D2::DeltaComplex; fix_points::Bool = true) :: Bool
    nv(D1) == nv(D2) && ne(D1) == ne(D2) && np(D1) == np(D2) || return false
    D = deepcopy(D1)
    fgv = FGVertex(D, fix_points)
    return is_isomorphic(fgv, D2, fix_points=fix_points)
end


"""
    is_isomorphic(D::DeltaComplex, D2::DeltaComplex; kwargs...) -> Bool

Return `true` if `D` is identical to `D2` up to a renaming of the edges and triFaces (and points if fix_points=false).

`D` is supposed to be already renamed in a canonical way. 

# Arguments
- `fix_points::Bool=true` : If is set to `false`, then the isomorphism would also allow a renaming of the points. 
"""
function is_isomorphic(fgv::FGVertex, D2::DeltaComplex; fix_points::Bool = true) :: Bool
    D = fgv.D
    numV = nv(D); numE = ne(D); numP = np(D)
    if numV != nv(D2) || numE != ne(D2) || numP != np(D2) || sort(point_degrees(D)) != sort(point_degrees(D2))
        return false
    end
    #point mapping
    A_tri = adjacency_matrix_triangulation(D)
    A_tri_2 = adjacency_matrix_triangulation(D2)
    if !fix_points
        permutations_points = mcKay_points(D2)
        if length(permutations_points) != fgv.num_point_perms
            return false
        end
        i = 1
        while i <= length(permutations_points)
            sig = invert_permutation(permutations_points[i])
            if !all(A_tri.== A_tri_2[sig, sig])
                deleteat!(permutations_points, i)
            else
                i += 1
            end
        end
        if isempty(permutations_points)
            return false
        end
    else
        if !all(A_tri .== A_tri_2)
            return false
        end
        permutations_points = [collect(1:numP)]
    end

    for perm_points in permutations_points
        D2_c = (length(permutations_points)>1 ? deepcopy(D2) : D2)
        rename_points!(D2_c, perm_points)
        
        #Triface mapping
        A_delta = adjacency_matrix_deltaComplex(D)
        A_delta_2 = adjacency_matrix_deltaComplex(D2_c)
        permutations_trifaces = mcKay_vertices(D2_c)
        if length(permutations_trifaces) != fgv.num_vertex_perms
            continue
        end
        i = 1
        while i <= length(permutations_trifaces) #TODO I think this check is pointless as either they are all valid or none are
            perm_trifaces = invert_permutation(permutations_trifaces[i])
            if !all(A_delta .== A_delta_2[perm_trifaces, perm_trifaces])
                deleteat!(permutations_trifaces, i)
            else
                i += 1
            end
        end
        if isempty(permutations_trifaces)
            continue
        end

        for perm_trifaces in permutations_trifaces
            D2_cc = (length(permutations_trifaces)>1 ? deepcopy(D2_c) : D2_c)
            rename_vertices!(D2_cc, perm_trifaces)

            #DualEdge mapping
            permutations_edges = mcKay_edges(D2_cc)
            if length(permutations_edges) != fgv.num_edge_perms
                continue
            end
            i = 1
            while i <= length(permutations_edges)
                perm_edges = invert_permutation(permutations_edges[i])
                if !all(j -> is_similar(D2_cc.E[perm_edges[j]], D.E[j]), 1:numE)
                    deleteat!(permutations_edges, i)
                else
                    i += 1
                end
            end
            if isempty(permutations_edges)
                continue
            end
            for perm_edges in permutations_edges
                D2_ccc = (length(permutations_edges)>1 ? deepcopy(D2_cc) : D2_cc)
                rename_edges!(D2_ccc, perm_edges) 
                bo = true
                for n in 1:numV
                    if !is_equivalent(D.V[n],D2_ccc.V[n])
                        return false
                    end
                end
                return true
            end
        end
    end
    return false
end

function is_equivalent(T1::TriFace, T2::TriFace)
    offset = -1
    e1_ids = collect(edges_id(T1))
    e2_ids = collect(edges_id(T2))
    if all(e1_ids .== e2_ids)
        offset = 0
    elseif all(e1_ids .== e2_ids[[2,3,1]])
        offset = 1
    elseif all(e1_ids .== e2_ids[[3,1,2]])
        offset = 2
    else
        return false
    end
    return all(i -> T1.points[i] == T2.points[(i+offset-1)%3 + 1] && e1_ids[i] == e2_ids[(i+offset-1)%3 + 1], 1:3) 
end



"""
    diameter(G::FlipGraph)

Return the diameter of the portion of the flipgraph that has been constructed.
"""
function diameter(G::FlipGraph)
    return diameter(adjacency_matrix(G.adjList))
end

"""
    split(V::Vector{<:Integer}, degs::Vector{<:Integer})

Split `V` into partitions according to their degrees from smallest to biggest.
"""
function split(V::Vector{<:Integer}, degs::Vector{<:Integer}) 
    sV = Vector{Vector{Int}}()
    unique_degs = sort!(unique(degs))
    j = 1 
    k = length(V)
    while k > 0
        W = Vector{Int}()
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
    mcKay_points(D::DeltaComplex; only_one::Bool=false)::Vector{Vector{Int}}

Apply a version of McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the points which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that Point 1 becomes Point `p[1]`, Point 2 becomes Point `p[2]`,...
If `only_one=true`, the algorithm stops after finding one valid permutation.
"""
function mcKay_points(D::DeltaComplex; only_one::Bool=false)::Vector{Vector{Int}}
    A = multi_adjacency_matrix_triangulation(D)

    #replace partitions by partitions as long as there are 2 elements in the same partition ...
    #...that may be differentiated by their relative degrees to another partition
    function makeEquitable!(p::Vector{Vector{T}}) where T<:Integer
        i = 1; j = 1
        while i <= length(p)
            rDegs = relative_degrees(A, p[i], p[j])
            if !all(rDegs .== rDegs[1]) 
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
    end

    function split_by_adjacent_holes!(p::Vector{Vector{T}}) where T<:Integer
        #compares the two arrays in a lexicographic sense returns true if u <= v (lexicographically)
        function is_less(u,v)
            for i in eachindex(u)
                if u < v
                    return true
                elseif u > v
                    return false
                end
            end
            return true
        end
        #points = collect()#reduce(vcat, collect((length(p[i])>1 ? p[i] : [])  for i in 1:length(p)))
        holecounts = [zeros(Int, length(D.holes)) for i in 1:D.num_points.x]
        for d in edges(D)
            x, y = triangle_edge(D, d)
            for c in D.edge_crossings[id(d)]
                holecounts[x][c.hole_id] += 1
                holecounts[y][c.hole_id] += 1
            end
        end
        i = 1
        while i <= length(p)
            if length(p[i]) > 1
                bigger_than = zeros(Int, length(p[i]))
                for j in eachindex(p[i])
                    for k in 1:j-1
                        if is_less(holecounts[p[i][j]], holecounts[p[i][k]])
                            bigger_than[k] += 1
                        else
                            bigger_than[j] += 1
                        end
                    end
                end
                newVs = split(p[i], bigger_than)
                #replace the old partition by the new ones
                popat!(p,i)
                for V in newVs
                    insert!(p, i, V)
                    i += 1
                end
            else
                i+=1
            end
        end
    end

    n = np(D)
    p = split(collect(1:n), degrees(A))
    makeEquitable!(p)

    if length(p) == n #there is only one canonical permutation
        return Vector{Vector{Int}}([invert_permutation(reduce(vcat, p))])
    else
        #split_by_adjacent_holes!(p)
        makeEquitable!(p)
    end
    if length(p) == n #|| only_one#there is only one canonical permutation
        return Vector{Vector{Int}}([invert_permutation(reduce(vcat, p))])
    end

    while only_one     
        i = 1
        while length(p[i]) == 1 #i = index of first partition that is not trivial
            i += 1 
        end 
        V = popat!(p, i)
        insert!(p, i, [popat!(V, 1)])
        insert!(p, i+1, V)
        makeEquitable!(p)
        if length(p) == n
            return Vector{Vector{Int}}([invert_permutation(reduce(vcat, p))])
        end
    end
    
    #split the first partition that has more than 2 elements 
    queue = Vector{Vector{Vector{Int}}}([p])
    leafs = Vector{Vector{Vector{Int}}}()
    while !isempty(queue)
        p = popfirst!(queue)
        i = 1       
        while length(p[i]) == 1; i += 1 end #i = index of first partition that is not trivial
        for j in 1:length(p[i]) #replace the i-th partition by isolating its j-th element 
            pp = deepcopy(p)
            V = popat!(pp, i)
            insert!(pp, i, [popat!(V, j)])
            insert!(pp, i+1, V)
            makeEquitable!(pp)
            if length(pp) != n
                push!(queue, pp)
            else
                push!(leafs, pp)
            end
        end
    end

    return [invert_permutation(reduce(vcat, p)) for p in leafs]   #Permutations
end

"""
    mcKay_vertices(D::DeltaComplex)::Vector{Vector{Int}}

Apply a version of McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the triFaces which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that TriFace 1 becomes TriFace p[1], TriFace 2 becomes TriFace p[2],...\\
If `only_one=true`, the algorithm stops after finding one valid permutation.

"""
function mcKay_vertices(D::DeltaComplex; only_one::Bool=false) ::Vector{Vector{Int}}
    np1 = np(D)
    np2 = np1*np1
    A = adjacency_matrix_deltaComplex(D)

    #replace partitions by partitions as long as there are 2 elements in the same partition ...
    #...that may be differentiated by their relative degrees to another partition
    function makeEquitable!(p::Vector{Vector{T}}) where T<:Integer
        i = 1; j = 1
        while i <= length(p)
            rDegs = relative_degrees(A, p[i], p[j])
            if !all(rDegs.==rDegs[1]) 
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
    end

    #give each combination of three points a distinct value (1,1,1)<(1,1,2)=(1,2,1)=(2,1,1)<(1,1,3)...  and (1,2,3)<(1,3,2)
    function pointValue(T::TriFace)
        p1,p2,p3 = T.points.-1
        return min(p1*np2 + p2*np1 + p3, p2*np2 + p3*np1 + p1, p3*np2 + p1*np1 + p2)
    end

    n = nv(D)
    p = split(collect(1:n), collect(pointValue(T) for T in D.V))
    makeEquitable!(p)

    if length(p) == n #|| only_one #there is only one canonical permutation
        return Vector{Vector{Int}}([invert_permutation(reduce(vcat, p))])
    end
    
    while only_one     
        i = 1
        while length(p[i]) == 1 #i = index of first partition that is not trivial
            i += 1 
        end 
        V = popat!(p, i)
        insert!(p, i, [popat!(V, 1)])
        insert!(p, i+1, V)
        makeEquitable!(p)
        if length(p) == n
            return Vector{Vector{Int}}([invert_permutation(reduce(vcat, p))])
        end
    end
    
    #split the first partition that has more than 2 elements 
    queue = Vector{Vector{Vector{Int}}}([p])
    leafs = Vector{Vector{Vector{Int}}}()
    while !isempty(queue)
        if length(queue)>=100
            println("Oh no")
        end
        p = popfirst!(queue)
        i = 1
        while length(p[i]) == 1; i += 1 end #i = index of first partition that is not trivial
        for j in eachindex(p[i]) #replace the i-th partition by isolating its j-th element 
            pp = deepcopy(p)
            V = popat!(pp, i)
            insert!(pp, i, [V[j]])
            popat!(V, j)
            insert!(pp, i+1, V)
            makeEquitable!(pp)
            if length(pp) != n
                push!(queue, pp)
            else
                push!(leafs, pp)
            end
        end
    end

    return [invert_permutation(reduce(vcat, p)) for p in leafs]   #sigma_pi's
end


"""
    mcKay_edges(D::DeltaComplex; only_one::Bool=false)::Vector{Vector{Int}}

Apply McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the triFaces which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that DualEdge 1 becomes DualEdge p[1], DualEdge 2 becomes DualEdge p[2],...\\
If `only_one=true`, the algorithm stops after finding one valid permutation.
"""
function mcKay_edges(D::DeltaComplex; only_one::Bool=false)::Vector{Vector{Int}}
    b1 = max(np(D), nv(D))
    b2 = b1*b1

    #give each combination of two vertices and 2 points a distinct value 
    function edgeValue(d::DualEdge)
        t1,t2 = d.triangles
        p1,p2 = points(D, d)
        return min(t1+t2*b1, t2+t1*b1)*b2 + min(p1+p2*b1, p2+p1*b1)
    end

    m = ne(D)
    p = split(collect(1:m), collect(edgeValue(d) for d in D.E))
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
        return Vector{Vector{Int}}([invert_permutation(reduce(vcat, p))])
    end
    
    #split the first partition that has more than 2 elements 
    queue = Vector{Vector{Vector{Int}}}([p])
    leafs = Vector{Vector{Vector{Int}}}()
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
                push!(leafs, pp)
            end
        end
    end

    return [invert_permutation(reduce(vcat, p)) for p in leafs]   #sigma_pi's
end
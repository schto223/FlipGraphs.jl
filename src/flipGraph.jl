export FlipGraph, construct_FlipGraph
export edges, has_edge, ne, nv, vertices
export diameter
export mcKay_points, mcKay_triFaces, mcKay_dualEdges
export is_isomorph, is_isomorph_to
export invert_perm

"""
    struct FlipGraph <: AbstractGraph{Int}

A Graph representing the FlipGraph of a `deltaComplex`.\\
Vertices are different triangulations of the same surface.\\
Two vertices are linked by an edge, if the respective graphs differ only by a single flip.
"""
struct FlipGraph <: AbstractGraph{Int}    
    V::Vector{HoleyDeltaComplex}
    adjList::Vector{Vector{Int}}

    function FlipGraph()
        new(Vector{HoleyDeltaComplex}(), Vector{Vector{Int}}())
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
    E = collect(Edge(i,j) for i = 1:length(G.V) for j in G.adjList[i])
    return filter!(e->(src(e)>dst(e)), E)
end 

edgetype(G::FlipGraph) = SimpleEdge{Int}
has_edge(G::FlipGraph, e::Edge) = (dst(e) ∈ G.adjList[src(e)])
has_edge(G::FlipGraph, s::Integer, d::Integer) = (d ∈ G.adjList[s])
has_edge(G::FlipGraph, HD1::HoleyDeltaComplex, HD2::HoleyDeltaComplex) = has_edge(G, findfirst(x->x==HD1, G.V), findfirst(x->x==HD2, G.V))
has_vertex(G::FlipGraph, v) = (1 <= v <= nv(G))
inneighbors(G::FlipGraph, v) = G.adjList[v]
ne(G::FlipGraph) = sum(size(G.adjList[i],1) for i=1:length(G.adjList))÷2
nv(G::FlipGraph) = length(G.V)
outneighbors(G::FlipGraph,v) = G.adjList[v]
vertices(G::FlipGraph) = G.V
is_directed(G::FlipGraph) = false
is_directed(::Type{FlipGraph}) = false


function add_edge!(G::FlipGraph, v, w) 
    if !has_edge(G, v, w) && v!=w
        push!(G.adjList[v],w)
        push!(G.adjList[w],v)
    end
end

function add_vertex!(G::FlipGraph, HD::HoleyDeltaComplex) 
    push!(G.V, HD)
    push!(G.adjList,[])
end

function remove_edge!(G::FlipGraph, e::Edge)
    deleteat!(G.adjList[src(e)], findfirst(x->x==dst(e), G.adjList[src(e)]))
    deleteat!(G.adjList[dst(e)], findfirst(x->x==src(e), G.adjList[dst(e)]))
end


"""
    construct_FlipGraph(HD::HoleyDeltaComplex, depth::Integer, modular::Bool=true)
    
Construct the **FlipGraph** for the HoleyDeltaComplex `HD`.  

If `modular` is true, then vertices of the FlipGraph are the classes of isomorphisms up to renaming the vertices. Each class is represented by one of its elements.\\
If `modular` is false, then each vertex is a different triangulation of the initial graph `g`.\\
By default, `modular` is set to `true`
"""
function construct_FlipGraph(HD::HoleyDeltaComplex, depth::Integer, fix_points::Bool = true)
    G = FlipGraph()
    if !fix_points
        HD = rename_points!(HD, mcKay_points(HD)[1])
    end
    add_vertex!(G, HD)
    queue = Vector{Tuple{HoleyDeltaComplex, Int, Int}}()  #(HD, index, depth)
    push!(queue, (HD,1,0))
    while !isempty(queue)
        HD, ind_HD, d = popfirst!(queue)
        for e in 1:ne(HD), left in (true) #TODO check if both flip directions are needed
            if is_flippable(HD, e)
                new_HD = flip(HD, e, left)
                is_new = true
                for i in eachindex(G.V)
                    if is_isomorph_to(G.V[i], new_HD, fix_points)
                        add_edge!(G, ind_HD, i)
                        is_new = false
                        break
                    end
                end
                if is_new #add the new DeltaComplex to the Graph
                    fix_points || rename_points!(new_HD, mcKay_points(new_HD)[1])
                    rename_vertices!(new_HD, mcKay_triFaces(new_HD)[1])
                    rename_edges!(new_HD, mcKay_dualEdges(new_HD)[1])
                    remove_holeloops!(new_HD)
                    add_vertex!(G, new_HD)
                    add_edge!(G, ind_HD, nv(G))
                    if d < depth
                        push!(queue, (G.V[end], nv(G), d+1))
                    end
                end
            end
        end
    end
    return G
end


"""
    is_isomorph(HD1::HoleyDeltaComplex, HD2::HoleyDeltaComplex, fix_points::Bool = true)

Return `true` if HD1 is isomorph to HD2 up to a renaming of the vertices, edges and if !fix_points also points.
"""
function is_isomorph(HD1::HoleyDeltaComplex, HD2::HoleyDeltaComplex, fix_points::Bool = true) :: Bool
    nv(HD1) == nv(HD2) && ne(HD1) == ne(HD2) && np(HD1) == np(HD2) || return false
    HD = deepcopy(HD1)
    remove_holeloops!(HD)
    fix_points || rename_points!(HD, mcKay_points(HD)[1])
    rename_vertices!(HD, mcKay_triFaces(HD)[1])
    rename_edges!(HD, mcKay_dualEdges(HD)[1])
    return is_isomorph_to(HD, HD2, fix_points)
end


"""
    is_isomorph_to(HD::HoleyDeltaComplex, HD2::HoleyDeltaComplex, fix_points::Bool = true)

Return true if `HD` is identical to `HD2` up to a renaming of the points, edges and triFaces.

`HD` is supposed to be already renamed in a canonical way. 

If `fix_points`=true, points are considered to be fixed and unchangeable. 

See also [`is_isomorph`](@ref)
"""
function is_isomorph_to(HD::HoleyDeltaComplex, HD2::HoleyDeltaComplex, fix_points::Bool = true)
    numV = nv(HD); numE = ne(HD); numP = np(HD)
    if numV != nv(HD2) || numE != ne(HD2) || numP != np(HD2) ||
        !fix_points && sort(point_degrees(HD)) != sort(point_degrees(HD2))
        return false
    end
    #point mapping
    A_tri = adjacency_matrix_triangulation(HD.D)
    A_tri_2 = adjacency_matrix_triangulation(HD2.D)
    if !fix_points
        permutations_points = mcKay_points(HD2)
        i = 1
        while i <= length(permutations_points)
            sig = invert_perm(permutations_points[i])
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
    remove_holeloops!(HD2)

    for perm_points in permutations_points
        HD2_c = deepcopy(HD2)
        rename_points!(HD2_c, perm_points)
        
        #Triface mapping
        A_delta = adjacency_matrix_deltaComplex(HD)
        A_delta_2 = adjacency_matrix_deltaComplex(HD2_c)
        permutations_trifaces = mcKay_triFaces(HD2_c)
        i = 1
        while i <= length(permutations_trifaces) #TODO I think this check is pointless as either they are all valid or none are
            perm_trifaces = invert_perm(permutations_trifaces[i])
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
            HD2_cc = deepcopy(HD2_c)
            rename_vertices!(HD2_cc, perm_trifaces)

            #DualEdge mapping
            permutations_edges = mcKay_dualEdges(HD2_cc)
            i = 1
            while i <= length(permutations_edges)
                perm_edges = invert_perm(permutations_edges[i])
                if !all(j -> length(HD2_cc.edge_crossings[perm_edges[j]]) == length(HD.edge_crossings[j]), 1:numE) ||
                       !all(j -> is_similar(HD2_cc.D.E[perm_edges[j]], HD.D.E[j]), 1:numE)      
                    deleteat!(permutations_edges, i)
                else
                    i += 1
                end
            end
            if isempty(permutations_edges)
                continue
            end
            #Check if there is a permutation for which the hole crossings are identical
            for perm_edges in permutations_edges
                HD2_ccc = deepcopy(HD2_cc)
                rename_edges!(HD2_ccc, perm_edges) 
                bo = true
                i = 1
                #Check if the crossings along the holes are identical
                while i <= length(HD.holes) && bo
                    C = get_crossing(HD.holes[i])
                    C2 = get_crossing(HD2_ccc.holes[i])
                    C2_start = C2
                    skipfirstcheck = true
                    bo = false
                    while (C2.key_holeposition != C2_start.key_holeposition || skipfirstcheck) && !bo
                        skipfirstcheck = false
                        bo = C2.edge_id == C.edge_id
                        if bo #found a valid starting position to align the holes
                            C2_next = C2.next
                            C_next = C.next
                            #check if all the following crossings align as well
                            while C2_next.key_holeposition != C2_start.key_holeposition && bo 
                                bo = C2_next.edge_id == C_next.edge_id
                                C2_next = C2_next.next
                                C_next = C_next.next
                            end
                        end
                        if !bo #C2 cannot be aligned with C
                            C2 = C2.next
                        end
                    end
                    i += 1
                end
                #check if the edge crossings match
                if bo 
                    i = 1
                    while i <= length(HD.edge_crossings) && bo
                        ec_holes = [c.hole_id for c in HD.edge_crossings[i]]
                        ec_holes_2 = [c.hole_id for c in HD2_ccc.edge_crossings[i]]
                        ec_directions = [c.going_in for c in HD.edge_crossings[i]]
                        ec_directions_2 = [c.going_in for c in HD.edge_crossings[i]]
                        if !(ec_holes == ec_holes_2 && ec_directions == ec_directions_2 ||
                            ec_holes == reverse(ec_holes_2) && all(ec_directions .!== reverse(ec_directions_2)))
                            bo = false
                        end
                        i += 1
                    end
                end
                if bo 
                    return true
                end
            end
        end
    end
    return false
end

"""
    invert_perm(p::Vector{Int})

Invert a permutation `p` such that p[i] = j becomes p[j] = i
"""
function invert_perm(p::Vector{Int})
    pp = zeros(Int, length(p))
    foreach( i -> pp[p[i]] = i ,eachindex(p))
    return pp
end

function diameter(G::FlipGraph)
    return diameter(adjacency_matrix(G.adjList))
end

"""
    split(V::Vector{<:Integer}, degs::Vector{<:Integer})

Split V into partitions according to their degrees from smallest to biggest
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
    mcKay_points(HD::HoleyDeltaComplex)::Vector{Vector{Int}}

Apply McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the points which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that Point 1 becomes Point p[1], Point 2 becomes Point p[2],...
"""
function mcKay_points(HD::HoleyDeltaComplex)::Vector{Vector{Int}}
    A = multi_adjacency_matrix_triangulation(HD.D)

    #replace partitions by partitions as long as there are 2 elements in the same partition ...
    #...that may be differentiated by their relative degrees to another partition
    function makeEquitable!(p::Vector{Vector{T}}) where T<:Integer
        i = 1; j = 1
        while i <= length(p)
            rDegs = relative_degrees(A, p[i], p[j])
            if !allequal(rDegs) 
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

    n = np(HD)
    p = split(collect(1:n), degrees(A)) #TODO replace point_degrees by function using adjacency matrix
    makeEquitable!(p)

    if length(p) == n #there is only one canonical permutation
        return Vector{Vector{Int}}([invert_perm(reduce(vcat, p))])
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

    return [invert_perm(reduce(vcat, p)) for p in leafs]   #Permutations
end

"""
    mcKay_triFaces(HD::HoleyDeltaComplex)::Vector{Vector{Int}}

Apply McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the triFaces which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that TriFace 1 becomes TriFace p[1], TriFace 2 becomes TriFace p[2],...
"""
function mcKay_triFaces(HD::HoleyDeltaComplex)::Vector{Vector{Int}}
    np1 = np(HD)
    np2 = np1*np1
    A = adjacency_matrix_deltaComplex(HD.D)

    #replace partitions by partitions as long as there are 2 elements in the same partition ...
    #...that may be differentiated by their relative degrees to another partition
    function makeEquitable!(p::Vector{Vector{T}}) where T<:Integer
        i = 1; j = 1
        while i <= length(p)
            rDegs = relative_degrees(A, p[i], p[j])
            if !allequal(rDegs) 
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

    n = nv(HD)
    p = split(collect(1:n), collect(pointValue(T) for T in HD.D.V))
    makeEquitable!(p)

    if length(p) == n #there is only one canonical permutation
        return Vector{Vector{Int}}([invert_perm(reduce(vcat, p))])
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

    return [invert_perm(reduce(vcat, p)) for p in leafs]   #sigma_pi's
end


"""
    mcKay_dualEdges(HD::HoleyDeltaComplex)::Vector{Vector{Int}}

Apply McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the triFaces which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that DualEdge 1 becomes DualEdge p[1], DualEdge 2 becomes DualEdge p[2],...
"""
function mcKay_dualEdges(HD::HoleyDeltaComplex)::Vector{Vector{Int}}
    b1 = max(np(HD), nv(HD))
    b2 = b1*b1

    #give each combination of two vertices and 2 points a distinct value 
    function edgeValue(d::DualEdge)
        t1,t2 = d.triangles
        p1,p2 = points(HD.D, d)
        return min(t1+t2*b1, t2+t1*b1)*b2 + min(p1+p2*b1, p2+p1*b1)
    end

    m = ne(HD)
    p = split(collect(1:m), collect(edgeValue(d) for d in HD.D.E))
    i = 1
    while i <= length(p)
        if length(p[i]) > 1
            pp = popat!(p,i)
            pps = split(pp, collect(length(HD.edge_crossings[pp[j]]) for j in eachindex(pp)))
            for ppp in pps
                insert!(p, i, ppp)
                i += 1
            end
        else
            i += 1
        end
    end
    #edges may only be in the same partition, if they share the same 2 endpoints and triangles
    #this is only possible for two edges at a time, they cannot be distinguished

    if length(p) == m #there is only one canonical permutation
        return Vector{Vector{Int}}([invert_perm(reduce(vcat, p))])
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
                push!(leafs, pp)
            end
        end
    end

    return [invert_perm(reduce(vcat, p)) for p in leafs]   #sigma_pi's
end
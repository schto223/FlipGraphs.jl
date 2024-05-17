export FlipGraph, construct_FlipGraph
export edges, has_edge, ne, nv, vertices
export diameter
export mcKay_points, mcKay_triFaces, mcKay_dualEdges
export is_isomorph, is_isomorph_to

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
        new(Vector{TriangulatedPolygon,1}(), Vector{Vector{Int}}())
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
has_edge(G::FlipGraph, s, d) = (d ∈ G.adjList[s])
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
    push!(G.V,g)
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
function construct_FlipGraph(HD::HoleyDeltaComplex, depth::Integer, fix_points::Bool=false)
    G = FlipGraph()
    fix_points || HD = rename_points!(HD, mcKay_points(HD)[1])
    add_vertex!(G, HD)
    
    queue = Vector{Tuple{HoleyDeltaComplex, Int, Int}}()  #(HD, index, depth)
    push!(queue, (HD,1,0))
    
    while !isempty(queue)
        HD, ind_HD, d = popfirst!(queue)
        for e in edges(HD), left in (true,false) #TODO check if both flip directions are needed
            if is_flippable(HD, e)
                new_HD = flip(HD, e, left)
                is_new = true
                for i = 1:length(G.V)
                    if is_isomorph_to!(G.V[i], new_HD, fix_points)
                        add_edge!(G, ind_HD, i)
                        is_new = false
                        break
                    end
                end
                if is_new #add the new DeltaComplex to the Graph
                    fix_points || rename_points!(new_HD, mcKay_points(new_HD)[1])
                    rename_trifaces!(new_HD, mcKay_triFaces(new_HD)[1])
                    rename_edges!(new_HD, mcKay_edges(new_HD)[1])
                    d += 1
                    add_vertex!(G, new_HD)
                    add_edge!(G, ind_HD, nv(G))
                    if d <= depth
                        push!(queue, (G.V[end], nv(G), d))
                    end
                end
            end
        end
    end
    return G
end


#"""
#    rename_vertices(g::TriGraph, sigma_pi::Array{Int,1})
#
#Rename the vertices of `g` by applying the permutation `sigma_pi`.
#"""
#function rename_vertices(g::TriangulatedPolygon, sigma_pi::Vector{<:Integer})
#    gg = TriangulatedPolygon(g.n, Vector{Vector{Int}}([[] for i in 1:g.n]))
#    for e in edges(g)
#        add_edge!(gg,sigma_pi[e.src], sigma_pi[e.dst])
#    end
#    return gg
#end

"""
    is_isomorph(HD1::HoleyDeltaComplex, HD2::HoleyDeltaComplex, fix_points::Bool = true)

Return `true` if HD1 is isomorph to HD2 up to a renaming of the vertices, edges and if !fix_points also points.
"""
function is_isomorph(HD1::HoleyDeltaComplex, HD2::HoleyDeltaComplex, fix_points::Bool = true) :: Bool
    nv(HD1) == nv(HD2) && ne(HD1) == ne(HD2) && np(HD1) == np(HD2) || return false
    HD = deepcopy(HD1)
    fix_points || rename_points!(HD, mcKay_points(HD)[1])
    rename_vertices!(HD, mcKay_triFaces(HD)[1])
    rename_edges!(HD, mcKay_edges(HD)[1])
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
    if !fix_points && sort(point_degrees(HD)) != sort(point_degrees(HD2)) || !fix_points && point_degrees(HD) != point_degrees(HD2)
        return false
    end
    #point mapping
    A_tri = adjacency_matrix_triangulation(HD.D)
    A_tri_2 = adjacency_matrix_triangulation(HD2.D)
    if !fix_points
        sigma_pi_points = mcKay_points(HD2)
        i = 1
        while i <= length(sigma_pi_points)
            sig = sigma_pi_points[i]
            if !all(A_tri.== A_tri_2[sig, sig])
                deleteat!(sigma_pi_points, i)
            else
                i += 1
            end
        end
        if isempty(sigma_pi_points)
            return false
        end
    else
        if !all(A_tri .== A_tri_2)
            return false
        end
        sigma_pi_points = collect(1:np(HD2))
    end


    
    for sigPi_points in sigma_pi_points
        HD2_c = deepcopy(HD2)
        rename_points!(HD2_c, sigPi_points)
        
        #Triface mapping
        A_delta = adjacency_matrix_deltaComplex(HD)
        A_delta_2 = adjacency_matrix_deltaComplex(HD2_c)
        sigma_pi_trifaces = mcKay_triFaces(HD2_c)
        i = 1
        while i <= length(sigma_pi_trifaces)
            sig = sigma_pi_trifaces[i]
            if !all(A_delta .== A_delta_2[sig, sig])
                deleteat!(sigma_pi_trifaces, i)
            else
                i += 1
            end
        end
        if isempty(sigma_pi_trifaces)
            return false
        end

        for sigPi_trifaces in sigma_pi_trifaces
            HD2_cc = deepcopy(HD2_c)
            rename_trifaces!(HD2_cc, sigPi_trifaces)

            #DualEdge mapping
            sigma_pi_edges = mcKay_dualEdges(HD2_cc)
            i = 1
            while i <= length(sigma_pi_edges)
                sigPi = sigma_pi_edges[i]
                if !all(HD2_cc.D.E[sigPi[j]] ≌ HD2.D.E[j] , j in 1:ne(HD)) #TODO define ≌(DualEdge,DualEdge)
                    deleteat!(sigma_pi_edges, i)
                else
                    i += 1
                end
            end
            if isempty(sigma_pi_edges)
                return false
            end

            for sigPi_edges in sigma_pi_edges
                HD2_ccc = deepcopy(HD2_cc)
                rename_edges!(HD2_ccc, sigPi_edges) 

                bo = true
                i = 1
                while i <= length(HD.holes) && bo
                    C = get_crossing(HD.holes[i])
                    C2 = get_crossing(HD_ccc.holes[i])
                    C2_start = C2
                    skipfirstcheck = true
                    bo = false
                    while (C2.key_holeposition != C2_start.key_holeposition || skipfirstcheck) && !bo
                        skipfirstcheck = false
                        bo = C2_start.edge_id == C.edge_id
                        if bo
                            C2_next = C2_start.next
                            C_next = C.next
                            while C2_next.key_holeposition != C2_start.key_holeposition && bo
                                bo = C2_next.edge_id == C_next.edge_id
                                C2_next = C2_next.next
                                C_next = C_next.next
                            end
                        end
                        if !bo
                            C2 = C2.next
                        end
                    end
                end
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
    deg = 0 
    k = length(V)
    while k > 0
        W = Vector{Int}()
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

"""
    mcKay_points(HD::HoleyDeltaComplex)::Vector{Vector{Int}}

Apply McKay's canonical graph labeling algorithm in order to determine all possible permutations 
of the points which give a canonical isomorphism class representant.

Return a vector of permutation vectors `p` such that Point 1 becomes Point p[1], Point 2 becomes Point p[2],...
"""
function mcKay_points(HD::HoleyDeltaComplex)::Vector{Vector{Int}}
    #replace partitions by partitions as long as there are 2 elements in the same partition ...
    #...that may be differentiated by their relative degrees to another partition
    function makeEquitable!(p::Vector{Vector{T}}, D::DeltaComplex) where T<:Integer
        i = 1; j = 1
        while i <= length(p)
            rDegs = relative_point_degrees(D, p[i], p[j])
            if !all(x -> x==rDegs[1], rDegs) 
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
    p = split(collect(1:n), point_degrees(HD.D))
    makeEquitable!(p, HD.D)

    if length(p) == n #there is only one canonical permutation
        return Vector{Vector{Int}}([reduce(vcat, p)])
    end
    
    #split the first partition that has more than 2 elements 
    queue = Vector{Vector{Vector{Int}}}([p])
    leafs = Vector{Vector{Vector{Int}}}()
    while !isempty(queue)
        p = popfirst!(queue)
        i = 1
        while length(p[i]) == 1; i += 1 end #i = index of first partition that is not trivial
        for j = 1:length(p[i]) #replace the i-th partition by isolating its j-th element 
            pp = deepcopy(p)
            V = popat!(pp, i)
            insert!(pp, i, [V[j]])
            popat!(V, j)
            insert!(pp, i+1, V)
            makeEquitable!(pp, HD.D)
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

    relative_TriFace_degrees(U::Vector{<:Integer}, V::Vector{<:Integer}) = [sum(A[u,V]) for u in U]
    
    #replace partitions by partitions as long as there are 2 elements in the same partition ...
    #...that may be differentiated by their relative degrees to another partition
    function makeEquitable!(p::Vector{Vector{T}}) where T<:Integer
        i = 1; j = 1
        while i <= length(p)
            rDegs = relative_TriFace_degrees(p[i], p[j])
            if !all(x -> x==rDegs[1], rDegs) 
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
        return Vector{Vector{Int}}([reduce(vcat, p)])
    end
    
    #split the first partition that has more than 2 elements 
    queue = Vector{Vector{Vector{Int}}}([p])
    leafs = Vector{Vector{Vector{Int}}}()
    while !isempty(queue)
        p = popfirst!(queue)
        i = 1
        while length(p[i]) == 1; i += 1 end #i = index of first partition that is not trivial
        for j = 1:length(p[i]) #replace the i-th partition by isolating its j-th element 
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

    n = nv(HD)
    p = split(collect(1:n), collect(edgeValue(d) for d in HD.D.E))
    #edges may only be in the same partition, if they share the same 2 endpoints and triangles
    #this is only possible for two edges at a time, they cannot be distinguished

    if length(p) == n #there is only one canonical permutation
        return Vector{Vector{Int}}([reduce(vcat, p)])
    end
    
    #split the first partition that has more than 2 elements 
    queue = Vector{Vector{Vector{Int}}}([p])
    leafs = Vector{Vector{Vector{Int}}}()
    while !isempty(queue)
        p = popfirst!(queue)
        i = 1
        while length(p[i]) == 1; i += 1 end #i = index of first partition that is not trivial
        for j = 1:length(p[i]) #replace the i-th partition by isolating its j-th element 
            pp = deepcopy(p)
            V = popat!(pp, i)
            insert!(pp, i, [V[j]])
            popat!(V, j)
            insert!(pp, i+1, V)
            if length(pp) != n
                push!(queue, pp)
            else
                push!(leafs, pp)
            end
        end
    end

    return [invert_perm(reduce(vcat, p)) for p in leafs]   #sigma_pi's
end
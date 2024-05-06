
export TriangulatedPolygon, flip!, mcKay, triangulatedPolygon 


"""
    struct TriangulatedPolygon <: AbstractGraph{Integer} 

A structure representing a triangulation of a convex polygon.
"""
struct TriangulatedPolygon <: AbstractGraph{Integer}    
    n::Int
    E::Array{Edge{Int}, 1}
    adjList::Array{Array{Int,1},1}

    function TriangulatedPolygon(n::Int)
        new(n, [], [[] for i=1:n])
    end

    function TriangulatedPolygon(n::Int, E::Array{Edge{Int}, 1}, adjList::Array{Array{Int,1},1})
        new(n,E,adjList)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", p::TriangulatedPolygon)
    print(io, string("TriangulatedPolygon with ",p.n, " vertices, and adjacency matrix:")); println(io)
    print(io, p.adjList)
end


edges(g::TriangulatedPolygon) = g.E
edgetype(g::TriangulatedPolygon) = SimpleEdge{Int}
has_edge(g::TriangulatedPolygon, e::Edge) = (e ∈ g.E)
has_edge(g::TriangulatedPolygon, s, d) = (d ∈ g.adjList[s])
has_vertex(g::TriangulatedPolygon, v) = (1 <= v && v <= g.n)
inneighbors(g::TriangulatedPolygon,v) = g.adjList[v]
ne(g::TriangulatedPolygon) = size(g.E, 1)
nv(g::TriangulatedPolygon) = g.n
outneighbors(g::TriangulatedPolygon,v) = g.adjList[v]
vertices(g::TriangulatedPolygon) = collect(1:g.n)
is_directed(g::TriangulatedPolygon) = false
is_directed(::Type{TriangulatedPolygon}) = false

function add_edge!(g::TriangulatedPolygon, v, w) 
    if !has_edge(g, v, w)
        push!(g.E, Edge(v,w))
        push!(g.adjList[v],w)
        push!(g.adjList[w],v)
    end
end

function remove_edge!(g::TriangulatedPolygon, e::Edge)
    deleteat!(g.E, findfirst(isequal(e), g.E))
    deleteat!(g.adjList[src(e)], findfirst(x->x==dst(e), g.adjList[src(e)]))
    deleteat!(g.adjList[dst(e)], findfirst(x->x==src(e), g.adjList[dst(e)]))
end


"""
    triangulatedPolygon(n::Int)

Create a triangulated convex n-gon.
"""
function triangulatedPolygon(n::Int)
    g = TriangulatedPolygon(n)

    for i = 1:n-1
        add_edge!(g, i, i+1)
    end
    add_edge!(g,n,1)
     
    if n <= 3
        return g
    end

    i = 2; j = n
    boo = true
    while i+1 < j
        add_edge!(g, i, j)
        if boo
            i+=1
        else
            j-=1
        end
        boo = !boo
    end

    return g
end

"""
    flip!(g::TriangulatedPolygon, e::Edge)

Flip the edge `e` in `g`.
"""
function flip!(g::TriangulatedPolygon, e::Edge)
    neigh1 = outneighbors(g, src(e))
    neigh2 = outneighbors(g, dst(e))
    S = intersect(neigh1, neigh2)
    u,v = S
    remove_edge!(g,e)
    add_edge!(g,u,v)
    return g
end

"""
    flip!(g::TriangulatedPolygon, e::Edge)

Return the triangulated polygon obtained by flipping the edge `e` in `g`.
"""
function flip(g::TriangulatedPolygon, e::Edge)  
    flip!(deepcopy(g), e) 
end

"""
    isFlippable(g::TriangulatedPolygon, e::Edge)

Return `true` if the edge `e` is flippable.\\

Note that for a triangulation of a convex polygon inner edges are always flippable, while outer edges cannot be flipped.
"""
function isFlippable(g::TriangulatedPolygon, e::Edge)
    neigh1 = outneighbors(g, e.src)
    neigh2 = outneighbors(g, e.dst)
    S = intersect(neigh1, neigh2)
    if size(S,1) < 2
        return false
    else
        return true
    end
end

function degrees(g::TriangulatedPolygon)
    return [size(g.adjList[i],1) for i=1:g.n]
end

len(v) = size(v,1)


mutable struct mcKayTreeNode #was mutable
    u::Int #::Base.RefValue{Int} #was just Int
    p::Array{Array{Int,1},1}
    children::Array{mcKayTreeNode}
    parent::mcKayTreeNode

    function mcKayTreeNode(n::Int)
        new(0,[collect(1:n)],[],nothing)
    end
    function mcKayTreeNode(u::Int,p::Array{Array{Int,1},1}, parent::mcKayTreeNode)
        new(u, p, [], parent)
    end
    function mcKayTreeNode(u::Int,p::Array{Array{Int,1},1})
        new(u, p, [] )
    end
end

function setParent!(child::mcKayTreeNode, parent::mcKayTreeNode)
    child.parent = parent
end

function addChild!(parent::mcKayTreeNode, child::mcKayTreeNode)
    push!(parent.children, child)
end


function relDegs(g::TriangulatedPolygon, U::Array{Int,1}, V::Array{Int,1})
    rdegs = zeros(Int, len(U))
    for i in 1:len(U), j in V
        if has_edge(g, U[i],j)
            rdegs[i]+=1
        end
    end
    return rdegs
end


#Apply the mcKay Algorithm in order to determine all possible permutations of the vertices which give a canonical isomorphism class representant.
function mcKay(g::TriangulatedPolygon)
    function split(V::Array{Int, 1}, degs::Array{Int, 1}) 
        sV = Array{Array{Int,1},1}()
        deg = 0 
        k = len(V)
        while k > 0
            W = Array{Int,1}()
            for i in 1:len(degs)
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

    function makeEquitable!(tNode::mcKayTreeNode, g::TriangulatedPolygon)
        p = tNode.p
        i = 1; j = 1
        while i <= len(p)
            rDegs = relDegs(g, p[i], p[j])
            if !all(x->x==rDegs[1],rDegs) 
                newVs = split(p[i], rDegs)
                popat!(p,i)
                for V in reverse(newVs)
                    insert!(p, i, V)
                end
                j = i
                i = 1
            else 
                j += 1
                if j > len(p)
                    j = 1
                    i += 1
                end
            end
        end    
    end

    p = split(collect(1:g.n), degrees(g))
    root = mcKayTreeNode(0, p)
    makeEquitable!(root, g)

    if len(root.p) == g.n
        return Array{Array{Int,1},1}([σ(root.p)])
    end

    queue = [root]::Array{mcKayTreeNode,1}
    leafs = Array{mcKayTreeNode,1}()
    while !isempty(queue)
        tNode = popfirst!(queue)::mcKayTreeNode
        i = 1
        while len(tNode.p[i]) == 1
            i += 1
        end
        for j = 1:len(tNode.p[i])
            child = deepcopy(tNode)
            addChild!(tNode, child)
            setParent!(child, tNode)
            empty!(child.children)
            V = popat!(child.p, i)
            insert!(child.p, i, [V[j]])
            popat!(V, j)
            insert!(child.p, i+1, V)
            makeEquitable!(child, g)
            if len(child.p)!=g.n
                push!(queue, child)
            else
                push!(leafs, child)
            end
        end
    end

    sigpis = Array{Array{Int,1},1}()
    for leaf in leafs
        push!(sigpis, σ(leaf.p))
    end

    return sigpis
            
end

# renamed from sigma
function σ(p::Array{Array{Int,1}, 1})
    sigpi = zeros(Int, len(p))
    for i = 1:len(p)
        sigpi[p[i][1]] = i
    end
    return sigpi
end




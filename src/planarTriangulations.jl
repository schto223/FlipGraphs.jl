
export TriGraph, flip!, mcKay, triangulatedPolygon 

struct TriGraph <: AbstractGraph{Integer}    
    n::Int
    E::Array{Edge{Int}, 1}
    adjList::Array{Array{Int,1},1}

    function TriGraph(n::Int)
        new(n, [], [[] for i=1:n])
    end

    function TriGraph(n::Int, E::Array{Edge{Int}, 1}, adjList::Array{Array{Int,1},1})
        new(n,E,adjList)
    end
end


edges(g::TriGraph) = g.E
edgetype(g::TriGraph) = SimpleEdge{Int}
has_edge(g::TriGraph, e::Edge) = (e ∈ g.E)
has_edge(g::TriGraph, s, d) = (d ∈ g.adjList[s])
has_vertex(g::TriGraph, v) = (1 <= v && v <= g.n)
inneighbors(g::TriGraph,v) = g.adjList[v]
ne(g::TriGraph) = size(g.E,1)
nv(g::TriGraph) = g.n
outneighbors(g::TriGraph,v) = g.adjList[v]
vertices(g::TriGraph) = collect(1:g.n)
is_directed(g::TriGraph) = false
is_directed(::Type{TriGraph}) = false

function add_edge!(g::TriGraph, v, w) 
    push!(g.E, Edge(v,w))
    push!(g.adjList[v],w)
    push!(g.adjList[w],v)
end

function remove_edge!(g::TriGraph, e::Edge)
    deleteat!(g.E, findfirst(isequal(e), g.E))
    deleteat!(g.adjList[src(e)], findfirst(x->x==dst(e), g.adjList[src(e)]))
    deleteat!(g.adjList[dst(e)], findfirst(x->x==src(e), g.adjList[dst(e)]))
end

function triangulatedPolygon(n::Int)
    g = TriGraph(n)

    for i = 1:n-1
        add_edge!(g, i, i+1)
    end
    add_edge!(g,n,1)
     
    if n <= 3
        return g
    end

    i = 2
    j = n
    b = true
    while i+1 < j
        add_edge!(g,i,j)
        if b
            i+=1
        else
            j-=1
        end
        b = !b
    end

    return g
end

function flip!(g::TriGraph, e::Edge)
    neigh1 = outneighbors(g, src(e))
    neigh2 = outneighbors(g, dst(e))
    S = intersect(neigh1, neigh2)
    u,v = S
    remove_edge!(g,e)
    add_edge!(g,u,v)
    return g
end

function flip(g::TriGraph, e::Edge)  
    flip!(deepcopy(g), e) 
end

function flippable(g::TriGraph, e::Edge)
    neigh1 = outneighbors(g, e.src)
    neigh2 = outneighbors(g, e.dst)
    S = intersect(neigh1, neigh2)
    if size(S,1) < 2
        return false
    else
        return true
    end
end


function degrees(g::TriGraph)
    return [size(g.adjList[i],1) for i=1:g.n]
end


len(v) = size(v,1)


function drawPNG(g::TriGraph, fName::String ="triGraph" )
    n = g.n
    a = 2*π/n 
    x = [cos(π/2 -a + a*i) for i = 1:n]
    y = [-sin(π/2 -a + a*i) for i = 1:n]
    nodeLabel = 1:n
    draw(PNG("img/"*fName * ".png", 1000px, 1000px), gplot(g, x, y, nodelabel=nodeLabel))
end

function plot(g::TriGraph)
    a = 2*π/n 
    x = [cos(π/2 -a + a*i) for i = 1:n]
    y = [-sin(π/2 -a + a*i) for i = 1:n]
    nodeLabel = 1:n
    gplot(g, x,y, nodelabel=nodeLabel)
end



mutable struct mcKayTreeNode
    u::Int
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


function relDegs(g::TriGraph, U::Array{Int,1}, V::Array{Int,1})
    rdegs = zeros(Int, len(U))
    for i in 1:len(U), j in V
        if has_edge(g, U[i],j)
            rdegs[i]+=1
        end
    end
    return rdegs
end

function mcKay(g::TriGraph)
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

    function makeEquitable!(tNode::mcKayTreeNode, g::TriGraph)
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
        return root
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
        push!(sigpis, sigma(leaf.p))
    end

    return sigpis
            
end


function sigma(p::Array{Array{Int,1}, 1})
    sigpi = zeros(Int, len(p))
    for i = 1:len(p)
        sigpi[p[i][1]] = i
    end
    return sigpi
end




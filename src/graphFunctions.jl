#import Graphs: adjacency_matrix
#using LinearAlgebra
export diameter, adjacency_matrix

function diameter(adjacency_matrix :: Matrix{<:Integer})
    n = size(adjacency_matrix,1)
    function Seidel(G::Matrix{<:Integer})
        if all(G[i,j]==1 || i==j for i in 1:n, j in 1:n)
            return G
        end
        
        H = Int.(G + G*G.>0) #- Matrix(I,n,n)
        foreach(i-> H[i,i]=0, 1:n)
        Dist = Seidel(H)
        #DG = Dist*G
        degrees = [sum(G[i,:]) for i=1:n]
        d = 2*Dist - (Dist*G .< Dist.*transpose(degrees))
        return d
    end

    d = Seidel(adjacency_matrix)
    return maximum(d)
end


function adjacency_matrix(adjList::Vector{Vector{T}}) :: Matrix{<:Integer} where {T<:Integer}
    n = size(adjList,1)
    A = zeros(Int,n,n)
    for i = 1:n
        for j in adjList[i]
            A[i,j]=1
        end
    end
    return A
    #rows = vcat([fill(i,size(adjList[i],1)) for i in 1:n]...)
    #cols = vcat(adjList...)
    #vals = ones(Int, size(cols,1))
    #return spzeros(rows, cols, vals)
end
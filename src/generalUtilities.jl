"""
    diameter(adjacency_matrix :: Matrix{<:Integer}) :: Int

Compute the diameter of a graph from its simple adjacency matrix.

All values in `adjacency_matrix` should be either 0 or 1.
"""
function diameter(adjacency_matrix :: Matrix{T}) :: T where T<:Integer
    n = size(adjacency_matrix, 1)
    function Seidel(G::Matrix{T})
        if all(G[i,j]==1 || i==j for i in 1:n, j in 1:n)
            return G
        end
        
        H = T.(G + G*G .> 0) #- Matrix(I,n,n)
        foreach(i -> H[i,i] = 0, 1:n)
        Dist = Seidel(H)
        degrees = [reduce(+, G[i,:]) for i in 1:n]
        d = 2*Dist - (Dist*G .< Dist.*transpose(degrees))
        return d
    end
    return maximum(Seidel(adjacency_matrix))
end

"""
    distances(adjacency_matrix :: Matrix{<:Integer}) :: Matrix{<:Integer}

Compute the shortest distance from any vertex to any other vertex in the graph for the given `adjacency_matrix`.

Return a `Matrix` whose entry at `(i,j)` is the length of a shortest path from `i` to `j`.
"""
function distances(adjacency_matrix :: Matrix{T}) :: Matrix{T} where T<:Integer
    n = size(adjacency_matrix, 1)
    function Seidel(G::Matrix{T}) :: Matrix{T}
        if all(G[i,j]==1 || i==j for i in 1:n, j in 1:n)
            return G
        end
        H :: Matrix{T} = T.(G + G*G .> 0) #- Matrix(I,n,n)
        foreach(i -> H[i,i] = 0, 1:n)
        Dist = Seidel(H)
        degrees = [reduce(+, G[i,:]) for i in 1:n]
        d :: Matrix{T} = 2*Dist - (Dist*G .< Dist.*transpose(degrees))
        return d
    end
    return Seidel(adjacency_matrix)
end

"""
    adjacency_matrix(adjList::Vector{Vector{<:Integer}}) :: Matrix{Int}

Construct the adjacency matrix from an adjacency list.
"""
function adjacency_matrix(adjList::Vector{Vector{T}}) :: Matrix{T} where {T<:Integer}
    n = size(adjList,1)
    A = zeros(T,n,n)
    for i = 1:n
        for j in adjList[i]
            A[i,j] = 1
        end
    end
    return A
end

"""
    invert_permutation(p::Vector{<:Integer})

Return the inverse of the permutation `p`. 

# Example
```julia-repl
julia> p = [2,1,4,5,3];
julia> p_inv = invert_perm(p); 
julia> show(p_inv)
[4, 2, 5, 1, 6, 3]
julia> show(p_inv[p])
[1, 2, 3, 4, 5, 6]
```
"""
function invert_permutation(p::Vector{<:T}) :: Vector{T} where T<:Integer
    pp = zeros(T, length(p))
    foreach( i -> pp[p[i]] = i ,eachindex(p))
    return pp
end


"""
    degrees(A::Matrix{<:Integer}) :: Vector{<:Integer}

Return a vector containing the degrees of every vertex given an adjacency matrix `A`.
"""
function degrees(A::Matrix{T}) ::Vector{T} where T<:Integer
    return reshape(sum(A, dims=2), size(A,1))
end


"""
    matrix_equal(A::Matrix{Int}, B::Matrix{Int}, p::Vector{Int}) :: Bool

returns `true` if `A == B[p,p]``
"""
function matrix_equal(A::Matrix{T}, B::Matrix{T}, p::Vector{<:Integer}) :: Bool  where T<:Integer     
    for i in size(A,1)
        for j in size(A,2)
            if B[p[i],p[j]] != A[i,j]
                return false
            end
        end
    end
    return true
end    

"""
    matrix_equal(A::Matrix{Int}, B::Matrix{Int}) :: Bool

Return `true` if `A` equals `B`. 

This function is much faster than calling `A==B`. However, `A` and `B` are assumed to have the same dimensions.
"""
function matrix_equal(A::Matrix{T}, B::Matrix{T}) ::Bool  where T<:Integer 
    for i in size(A,1)
        for j in size(A,2)
            if B[i,j] != A[i,j]
                return false
            end
        end
    end
    return true
end

"""
    relative_degrees(A::Matrix{<:Integer}, U::Vector{<:Integer}, V::Vector{<:Integer}) :: Vector{Int32}

Compute the relative degrees of points in `U` onto the subset of points `V`.
"""
function relative_degrees(A::Matrix{T}, U::Vector{<:Integer}, V::Vector{<:Integer}) :: Vector{T} where T<:Integer
    rel_degs = zeros(T, length(U))
    for i in eachindex(U)
        for vj in V
            rel_degs[i] += A[U[i], vj]
        end
    end
    return rel_degs
end

"""
    relative_degree(A::Matrix{<:Integer}, u::Integer, V::Vector{<:Integer}) :: Int32

Compute the number of edges going from `u` into any vertex in the subset of points `V`.

# Arguments 
-`A::Matrix{<:Integer}`: the adjacency matrix. `A[i,j] = 1` if there is an edge going from `i` to `j`
"""
function relative_degree(A::Matrix{T}, u::Integer, V::Vector{<:Integer}) :: T where T<:Integer
    rel_deg :: T = T(0)
    for vj in V
        rel_deg += A[u, vj]
    end
    return rel_deg
end
# General utilities

These are some methods that you might want to use, but they aren't specific to one of the main structures in this package.

```@docs
diameter(::Matrix{T}) where T<:Integer
distances
adjacency_matrix(::Vector{Vector{T}}) where T<:Integer
matrix_equal
invert_permutation
degrees
relative_degree
relative_degrees(::Matrix{T}, ::Vector{<:Integer}, ::Vector{<:Integer}) where T<:Integer
```

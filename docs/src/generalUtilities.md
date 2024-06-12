# General utilities

These are some methods that you might want to use, but they aren't specific to one of the main structures in this package.

```@docs
diameter(::Matrix{<:Integer})
distances
adjacency_matrix(::Vector{Vector{T}}) where T<:Integer
matrix_equal
invert_permutation
degrees(A::Matrix{<:Integer})
relative_degree
relative_degrees(::Matrix{<:Integer}, ::Vector{<:Integer}, ::Vector{<:Integer})
```
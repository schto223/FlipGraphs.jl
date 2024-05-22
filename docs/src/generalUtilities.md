# General utilities

These are some helping methods which youmight want to use, but aren't specific to one of the main structures in this package.

```@docs

diameter(::Matrix{<:Integer})
adjacency_matrix(::Vector{Vector{T}}) where T<:Integer
invert_permutation(::Vector{<:Integer})
degrees(A::Matrix{<:Integer})
```
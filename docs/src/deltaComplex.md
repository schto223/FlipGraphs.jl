# DeltaComplexes

## Modelling a triangulation on a closed surface

A **Δ-complex** is a representation of a triangulation on a closed surface.
To define a triangulation on a closed surface, it does not suffice to take into account vertices and edges. We will also need to take into account the triangular faces between them. Therefore the triangulations are modeled using an extension of their dual graph.

Vertices are triangular faces which in turn consist of three points, and three edges. These points and edges are not necessarily distinct.
Edges in the dual Graph(i.e. the Δ-complex) connect two triangular faces if they in turn share an edge. In order to avoid confusion between the edges of the triangulation and the edges in the dual graph, I will, hence forwards call the latter *dual edge*.

```@docs
    DeltaComplex
    DualEdge
    TriFace
```

### Construction of a DeltaComplex

This Module comes with some handy and easy to use tools to construct a triangulation of a surface:

```@docs
    delta_complex
    subdivide!(::DeltaComplex, ::Integer)
```
### Common Graph-like methods

`DeltaComplex` is not an implementation of Graphs.AbstractGraph. <!--[`Graphs.AbstractGraph`](@extref).--> However, as it is similar to a simple graph, I used the same notation and function names for simplicity

```@docs 
    np(::DeltaComplex)
    nv(::DeltaComplex)
    ne(::DeltaComplex)
    vertices(::DeltaComplex)
    get_vertex(::DeltaComplex, ::Integer)
    get_vertex_id(::DualEdge, ::Integer)
    edges(::DeltaComplex)
    get_edge(::DeltaComplex, ::Integer)
    get_edge(::DeltaComplex, ::Integer, ::Integer)
    id(::DualEdge)

    points(::TriFace)
    get_point(::TriFace, ::Integer)
    points(::DeltaComplex, ::DualEdge)

    vertices(::DeltaComplex, ::DualEdge)
    vertices_id(::DualEdge)

    get_edge(::TriFace, ::Integer) 
    get_edge_id(::TriFace, ::Integer)
    edges(::TriFace)
    edges(::DeltaComplex, ::Integer)
    edges_id(::DeltaComplex, ::Integer) 
    id(::TriFace) 

    other_endpoint(::DualEdge, ::Integer, ::Integer)
    is_similar(::DualEdge, ::DualEdge)
```

## Classifying the triangulation

Here are some usefull methods, to pull out general information about the Δ-Complex, and the triangulation it represents:

```@docs
    euler_characteristic(::DeltaComplex)
    genus(::DeltaComplex)
    diameter_triangulation(::DeltaComplex)
    diameter_deltaComplex(::DeltaComplex)
    adjacency_matrix_deltaComplex(::DeltaComplex)
    adjacency_matrix_triangulation(::DeltaComplex)
    multi_adjacency_matrix_triangulation(::DeltaComplex)
```

```@docs
    point_degrees(::DeltaComplex)
    relative_point_degrees(::DeltaComplex, ::Vector{<:Integer}, ::Vector{<:Integer})
```

```@docs
    rename_edges!(::DeltaComplex, ::Vector{<:Integer})
    rename_points!(D::DeltaComplex, p::Vector{<:Integer})
    rename_vertices!(::DeltaComplex, ::Vector{<:Integer})
```

## Non-orientable closed surfaces

Regarding on which way edges are glued together, the resulting surface may be non-orientable. These surfaces can also be modeled by `DeltaComplex` and all of the methods above (with the exception of [`genus`](@ref)) may still be applied.

```@docs
    is_orientable(::DeltaComplex)
    demigenus
    twist_edges!
```

## Flipping

A **flip** is defined as the action of replacing an edge in the triangulation by the other diagonal of the quadrilateral formed by the two triangles adjacent to the edge.

It has been shown, that the flipgraph of any closed surface is connected, hence it is possible, to obtain any triangulation by a finite number of flips.

```@docs
    flip!(::DeltaComplex, ::Integer)
    flip!(::DeltaComplex, ::DualEdge; ::Bool)
    is_flippable(::DeltaComplex, ::Integer)
    is_flippable(::DualEdge)
    random_flips!(::DeltaComplex, ::Integer)
    randomize!(::DeltaComplex; ::Integer, ::Integer, ::Integer, ::Integer)
```
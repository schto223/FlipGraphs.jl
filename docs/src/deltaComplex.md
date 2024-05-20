# DeltaComplexes

## Modelling a triangulation on a closed surface

A Δ-complex is a representation of a triangulation on a closed surface.
To define a triangulation on a closed surface, it does not suffice to take into account vertices and edges. We will also need to take into account the triangular faces between them. Therefore the triangulations are modeled using an extension of their dual graph.

Vertices are triangular faces which in turn consist of three points, and three edges. These points and edges are not necessarily distinct.
Edges in the dual Graph(i.e. the Δ-complex) connect two triangular faces if they in turn share an edge. In order to avoid confusion between the edges of the triangulation and the edges in the dual graph, I will, hence forwards call the latter *dual edge*

```@docs
    DeltaComplex
    DualEdge
    TriFace
    delta_complex
    subdivide!(::DeltaComplex, ::Integer)
    rename_edges!(::DeltaComplex)
    rename_points!(::DeltaComplex)
    rename_vertices!(::DeltaComplex)
    np(::DeltaComplex)
    nv(::DeltaComplex)
    ne(::DeltaComplex)
    points(::DeltaComplex)
    get_point(::DeltaComplex)
    vertices(::DeltaComplex)
    get_vertex(::DeltaComplex)
    edges(::DeltaComplex)
    get_edge(::DeltaComplex)
    edges_id(::DeltaComplex)
    get_edge_id(::DeltaComplex)
    id(::DeltaComplex)
    other_endpoint(::DeltaComplex)
    _is_similar(::DeltaComplex)

    euler_characteristic(::DeltaComplex)
    genus(::DeltaComplex)
    diameter_triangulation(::DeltaComplex)
    diameter_deltaComplex(::DeltaComplex)
    adjacency_matrix_deltaComplex(::DeltaComplex)
    adjacency_matrix_triangulation(::DeltaComplex)
    multi_adjacency_matrix_triangulation(::DeltaComplex)
    point_degrees(::DeltaComplex)
    relative_point_degrees(::DeltaComplex)

    flip!(::DeltaComplex)
    is_flippable(::DeltaComplex)
    random_flips!(::DeltaComplex)
    randomize!(::DeltaComplex)
```

## Non-orientable closed surfaces

Regarding on which way edges are glued together, the resulting surface may be non-orientable. These surfaces can also be modeled by `DeltaComplex` and allof the methods above (with the exception of [`genus`](@ref)) may still be applied.

```@docs
    is_orientable(::DeltaComplex)
    demigenus
    twist_edges!
```
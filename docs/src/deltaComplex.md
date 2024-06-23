# Triangulations of Closed Surfaces
```@meta
DocTestSetup = quote
    using FlipGraphs
end
```
A triangulation of a closed surface is composed of points on the closed surface, which are connected by a maximal number of arcs (isotopy classes of curves on the surface that start and end at fixed points without containing any loops) in such a way, that no two arcs are homotopic to each other.\
As these become very complex objects that are difficult and computationally complex to model, it is often more useful to look at the dual of a triangulation. This dual is called a *Δ-complex*. 

## DeltaComplex

A **Δ-complex** is a representation of a triangulation on a closed surface.
To define a triangulation on a closed surface, it does not suffice to consider vertices and edges. We will also need to consider the triangular faces between them. Therefore, the triangulations are modeled using an extension of their dual graph.

Vertices are triangular faces, which in turn consist of three points and three edges. These points and edges are not necessarily distinct.
Edges in the dual (i.e. the Δ-complex) connect two triangular faces if they, in turn, share an edge. To avoid confusion between the edges of the triangulation and the edges of the dual, I will henceforth refer to the latter as the *dual edge*.

### Structures

```@docs
    DeltaComplex
    TriFace
    DualEdge
```

### Construction of a DeltaComplex

This Module comes with some handy and easy to use tools to construct a triangulation of a surface:

```@docs
    deltacomplex
    deltacomplex_non_orientable
    subdivide!(::DeltaComplex, ::Integer)
```
### Extracting information from `DeltaComplex`'

`DeltaComplex` is not an implementation of *Graphs.AbstractGraph*. However, as it is similar to a simple graph, the same notation and function names were used for simplicity.

```@docs 
    np(::DeltaComplex)
    nv(::DeltaComplex)
    ne(::DeltaComplex)
    vertices(::DeltaComplex)
    get_vertex(::DeltaComplex, ::Integer)
    
    edges(::DeltaComplex)
    get_edge(::DeltaComplex, ::Integer)
    get_edge(::DeltaComplex, ::Integer, ::Integer)

    quadrilateral_edges(::DeltaComplex, ::DualEdge)

    is_similar(::DualEdge, ::DualEdge)
```

### Extracting information from `TriFace`s

`TriFace`s form the vertices of a `DeltaComplex`. They are defined by their 3 edges. Additionally, they store the 3 points that form their corners.

As the `TriFace`s are stored in a list in `DeltaComplex`, each `TriFace` has its proper `id` which is the same as its index in the list.

```@docs
    id(::TriFace)
```

If you wish to get the edges of a `TriFace`, you might want to use any of the following methods:

```@docs 
    get_edge(::TriFace, ::Integer) 
    get_edge_id(::TriFace, ::Integer)
    edges(::TriFace)    
    edges_id(::TriFace)
    edges(::DeltaComplex, ::Integer)
    edges_id(::DeltaComplex, ::Integer) 
```

If you wish to get the points that form the corners of a `TriFace`, you might want to use any of the following methods:

```@docs 
    has_point(::TriFace, ::Integer)
    get_point(::TriFace, ::Integer)
    points(::TriFace)
    points(::DeltaComplex, ::DualEdge)
```

### Extracting information from `DualEdge`s

`DualEdge`s form the edges in a `DeltaComplex`. They are defined by the 2 triangles that form the endpoints and the respective sides through which they connect them. In addition, they may be twisted. (See [Non-orientable closed surfaces](@ref) for more on twisted edges)
Contrary to a normal graph, it is important on which side of a `TriFace` an edge is incident to.

Like the `TriFace`s, `DualEdge`s are also stored in a list in `DeltaComplex`, each `DualEdge` has its proper `id` which is the same as its index in the list.

```@docs     
    id(::DualEdge)
```

If you wish to get the `TriFace`s that border the `DualEdge`, you might want to use any of the following methods:

```@docs     
    vertices(::DeltaComplex, ::DualEdge)
    vertices_id(::DualEdge)    
    get_vertex_id(::DualEdge, ::Integer)
    sides 
    get_side(::DualEdge, ::Integer)   
    other_endpoint(::DualEdge, ::Integer, ::Integer)
```

## Classifying the triangulation

Here are some useful methods, to pull out general information about the Δ-Complex, and the triangulation it represents:

```@docs    
    genus(::DeltaComplex)
    euler_characteristic(::DeltaComplex)
    diameter_triangulation(::DeltaComplex)
    diameter_deltaComplex(::DeltaComplex)
    diameter(::DeltaComplex)
    adjacency_matrix_deltacomplex(::DeltaComplex)
    multi_adjacency_matrix_deltacomplex(::DeltaComplex)
    adjacency_matrix_triangulation(::DeltaComplex)
    multi_adjacency_matrix_triangulation(::DeltaComplex)
```

```@docs
    point_degrees(::DeltaComplex)
```

### Relabeling

If you want to relabel/reorder the points, vertices or edges, you may do so, by providing a permutation vector `` p=[p_1, p_2, \ldots , p_n] `` which will relabel ``i`` as ``p_i``. 

```@docs
    rename_edges!(::DeltaComplex, ::Vector{<:Integer})
    rename_points!(::DeltaComplex, ::Vector{<:Integer})
    rename_vertices!(::DeltaComplex, ::Vector{<:Integer})
```

## Non-orientable closed surfaces

Regarding on how different sides of triangles are associated to each other, the resulting surface may be non-orientable. These surfaces can also be modeled by `DeltaComplex`, and all the methods above (except for [`genus`](@ref)) may still be applied.

```@docs
    is_orientable(::DeltaComplex)
    demigenus
    twist_edges!
    is_twisted
```

## Flipping

A **flip** is defined as the action of replacing an edge in the triangulation by the other diagonal of the quadrilateral formed by the two triangles adjacent to the edge.

It has been shown, that the flip graph of any closed surface is connected. Hence, it is possible, to obtain any triangulation by a finite number of flips.

```@docs
    flip(::DeltaComplex, ::Integer)
    flip!(::DeltaComplex, ::Integer)
    flip!(::DeltaComplex, ::DualEdge; ::Bool)
    is_flippable(::DeltaComplex, ::Integer)
    is_flippable(::DualEdge)
    random_flips!(::DeltaComplex, ::Integer)
    randomize!(::DeltaComplex; ::Integer, ::Integer, ::Integer, ::Integer)
```
# HoleyDeltaComplex

!!! warning "Warning"
    This extension of a deltaComplex does not achieve its goal of allowing to accurately deduce if two triangulations are the same.
    Although it is more accurate than the simple deltaComplex, there are still many cases, where two triangulations are deemed to be identical when they are in fact different

The [`DeltaComplex`](@ref) structure is usefull for many applications. However, it has a slight issue that two different triangulations of the same surface, may result in the same DeltaComplex. In other words, a DeltaComplex misses some information about the underlying triangulation and can therefore not be used if one would, for example, like to construct the flipgraph (or rather a portion of latter) of the triangulations on a given surface.

The [`HoleyDeltaComplex`](@ref) structure solves this issue, by keeping track of the holes through which every single edge passes, and the orders thereof. This comes of course with a decrease in efficiency. It is therefore advised to use the simpler [`DeltaComplex`](@ref) if this bonus information is not needed. Another drawback is, that the `HoleyDeltaComplex` structure currently does not work for non-orientable surfaces.


## Structures

```@docs
    HoleyDeltaComplex
    Hole
    Crossing
```

## Construction
```@docs
    holey_delta_complex
    subdivide!(::HoleyDeltaComplex, ::Integer)
```

### Common Graph-like methods

As the `HoleyDeltaComplex` is essentially an extension of the `DeltaComplex` structure, all methods which take DeltaComplexes also work with HoleyDeltaComplexes.

```@docs
    np(::HoleyDeltaComplex)
    nv(::HoleyDeltaComplex)
    ne(::HoleyDeltaComplex)

    vertices(::HoleyDeltaComplex) 
    vertices(::HoleyDeltaComplex, ::DualEdge) 
    get_vertex(::HoleyDeltaComplex, ::Integer)

    edges(::HoleyDeltaComplex)
    edges(::HoleyDeltaComplex, ::Integer)
    get_edge(::HoleyDeltaComplex, ::Integer)
    get_edge(::HoleyDeltaComplex, ::Integer,::Integer)
    edges_id(::HoleyDeltaComplex, ::Integer)

    points(::HoleyDeltaComplex, ::DualEdge)

```    
```@docs
    rename_edges!(::HoleyDeltaComplex, ::Vector{<:Integer})
    rename_points!(::HoleyDeltaComplex, ::Vector{<:Integer})
    rename_vertices!(::HoleyDeltaComplex, ::Vector{<:Integer})
```


## Classifying the triangulation

```@docs
    euler_characteristic(::HoleyDeltaComplex)
    genus(::HoleyDeltaComplex)
    diameter_triangulation(::HoleyDeltaComplex)
    diameter_deltaComplex(::HoleyDeltaComplex)
    adjacency_matrix_deltaComplex(::HoleyDeltaComplex)
    adjacency_matrix_triangulation(::HoleyDeltaComplex)
    multi_adjacency_matrix_triangulation(::HoleyDeltaComplex)
    point_degrees(::HoleyDeltaComplex)
    relative_point_degrees(::HoleyDeltaComplex, ::Vector{<:Integer}, ::Vector{<:Integer})
```

## Flipping

```@docs
    flip(::HoleyDeltaComplex, ::Integer)
    flip(::HoleyDeltaComplex, ::DualEdge)
    flip!(::HoleyDeltaComplex, ::Integer)
    flip!(::HoleyDeltaComplex, ::DualEdge)
    is_flippable(::HoleyDeltaComplex, ::Integer)
    random_flips!(::HoleyDeltaComplex, ::Integer)
```

## HoleyDeltaComplex only

```@docs
    num_crossings(::Hole)
    edge_crossings(::HoleyDeltaComplex, ::Integer)
    get_crossing(::Hole)
    remove_holeloops!(::HoleyDeltaComplex)
```
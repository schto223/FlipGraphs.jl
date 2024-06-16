# Modular Flip Graphs of Closed Orientable Surfaces
```@meta
DocTestSetup = quote
    using FlipGraphs
end
```
Contrary to flip graphs of planar triangulations like that of a convex polygon, the flip graph of a closed surface is generally infinitely large.
Therefore, it is impossible to construct the whole flip graph. Not only are they infinitely large, but they also "grow" very fast. For this reason, it is more interesting and achievable to look at the **modular flip graphs**. These are graphs whose vertices are homotopy classes of triangulations on closed surfaces.

It turns out that the vertices of modular flip graphs are equivalent to the *Δ-coplexes* which we have already modeled in [DeltaComplex](@ref).

## Structures

A **`FlipGraph`** consists of vertices that represent homotopy classes of `DeltaComplex`es and edges between them. 
In order to be able to efficiently compute the whole flip graph, vertices are not `DeltaComplex`es but have their own structure (`FGVertex`) which consists of a `DeltaComplex` with additional information.

```@docs
    FlipGraph
    FGVertex
```

## Constructors

Currently, it is only possible to create the modular flip graph of a closed *orientable* surface. Non-orientable surfaces would need an adapted approach to solving, as the way they are modeled right now would require a lot more options to be checked just to determine if two `DeltaComplex`es are isomorphic to each other.

As even modular flip graphs become very large quite quickly, you have the possibility of only building a local portion of a flip graph. The way this works is, by giving it a `DeltaComplex`, which forms the root vertex, and setting the `depth`, which will limit the vertices to only include those that are at a distance less or equal to this `depth` from the root.

```@docs
    flipgraph_modular
```

## Graph methods

As with [`FlipGraphPlanar`](@ref), there are a bunch of methods which overload some of the main functions from the [Graphs.jl](https://juliagraphs.org/Graphs.jl/stable/) package.

```@docs
    nv(::FlipGraph)
    ne(::FlipGraph)
    vertices(::FlipGraph)
    vertices_deltacomplex
    edges(::FlipGraph)
    has_vertex(::FlipGraph, ::Integer)
    get_vertex(::FlipGraph, ::Integer)
    has_edge(::FlipGraph,::Edge)
    has_edge(::FlipGraph,::Integer, ::Integer)    
    has_edge(::FlipGraph,::FGVertex, ::FGVertex)
    neighbors(::FlipGraph, ::Integer)
    diameter(::FlipGraph)
```

## Comparing Triangulations

One problem in deciding if two triangulations are equivalent, is that the naming of the vertices, edges and points is completely arbitrary.
In the act of flipping, `DualEdge`s and `TriFace`s "move around". It is therefore possible to obtain two `DeltaComplex`es representing the same triangulation of a surface but with different labeled `TriFace`s and `DualEdge`s. 

The following method will determine if two `DeltaComplex`es are isomorphic to each other and therefore represent the same vertex in the modular flip graph:

```@docs
    is_isomorphic(::DeltaComplex, ::DeltaComplex)
```

## Canonical labeling of `DeltaComplex`es

Checking every possible permutation of `TriFace` and `DualEdge` labeling was not an option, as the number of possibilities would blow up immediately. 

What one can do instead is try to find a canonical labeling. This module uses a version of *McKay's canonical graph labeling algorithm*[^1] to try and determine a unique labeling based on the relationship of vertices and edges to other vertices and edges. In general, it is not always possible to determine a unique labeling. However, with this method, we can reduce the number of labelings to a manageable number.

The `mcKay` methods each return a permutation vector `p` which can be interpreted as a Cauchy's one-line notation for permutations. 
For example, `p = [3,5,1,2,6,4]` would correspond to the following permutation:

```math
σ = \begin{pmatrix}
1 & 2 & 3 & 4 & 5 & 6\\
3 & 5 & 1 & 2 & 6 & 4
\end{pmatrix} = (1\ 3)(2\ 5\ 6\ 4)
```

This algorithm has been adapted to find canonical labelings for points, vertices and edges. These are not independent from each other. In fact, only the point labelings depend on themselves. mcKay_vertices will give different results basedon the point labelings, and mcKay_edges will give different results based on both point and vertex labelings. Therefore, they are to be used in sequence.

```@docs
    mcKay_points
    mcKay_vertices
    mcKay_edges
```


[^1]: Hartke, S.G., & Radcliffe, A.J. (2008). McKay ’ s Canonical Graph Labeling Algorithm.
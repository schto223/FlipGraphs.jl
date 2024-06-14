
# Flip Graphs of Convex Polygons

**Flip Graphs**, are obtained by considering triangulations on a fixed set of points as vertices, and connecting two vertices, 
if it is possible to get from one triangulation to the other by flipping a single edge.

## Structures

The main structure for this part of the package is the **`FlipGraphPlanar`** structure. 
This implements the `AbstractGraph` interface from [Graphs.jl](https://juliagraphs.org/Graphs.jl/stable/). 
It is therefore possible to use it with other packages that work with *Graphs.jl*.

```@docs
    FlipGraphPlanar
```

For computational purposes, the vertices of a `FlipGraphPlanar` object are not `TriangulatedPolygon`s. They have their own proper structure, which is composed of a `TriangulatedPolygon` object together with some helpful information about the triangulation. \
In the case of a flip graph of triangulations with labeled points, the `TriangulatedPolygon` is *the* unique triangulation of its kind.\ 
If the flip graph is, however, of triangulations with unlabeled points, then the `TriangulatedPolygon` is only one of the triangulations where the points have been labeled by one of the possible canonical labelings.

```@docs
    FGPVertex
```

## Constructors

To construct the flip graph of a convex polygon, you can either start from a `TriangulatedPolygon` or just set the number of vertices of the triangulated convex polygon.

```@docs
    flipgraph(::TriangulatedPolygon)
    flipgraph_planar
```

## Graph methods

The following methods overload some of the main functions from the [Graphs.jl](https://juliagraphs.org/Graphs.jl/stable/) package.

```@docs
    nv(::FlipGraphPlanar)
    ne(::FlipGraphPlanar)
    vertices(::FlipGraphPlanar)
    edges(::FlipGraphPlanar)
    has_vertex(::FlipGraphPlanar,::Integer)
    has_vertex(::FlipGraphPlanar,::FGPVertex)
    has_edge(::FlipGraphPlanar,s,d)
    has_edge(::FlipGraphPlanar,::Edge)

    neighbors(::FlipGraphPlanar, ::Integer)
    diameter(::FlipGraphPlanar)
```

## Comparing Triangulations

In order to compute the *flip graph*, one needs to be able to determine if two triangulations are identical (if the points are labeled) or isomorphic to each other (if the points are unlabeled). \
The first case is fairly simple, as two triangulations are identical if their adjacency lists are. The only difficulty here is that the order in the list of neighbors is not fixed.\
The second case is more challenging, as there are $n!$ different ways to label $n$ points. The way it is done in this package is to use a variation of *McKay's canonical graph labeling algorithm*[^1] in order to rattle the number of possible labelings down to a relatively small number.

```@docs    
    is_isomorphic(::TriangulatedPolygon, ::TriangulatedPolygon)
    is_isomorphic(::FGPVertex, ::TriangulatedPolygon, ::Array{Vector{T},1}) where T<:Integer
```

The following methods are used to build the flip graph; However, they can also be useful elsewhere:

```@docs    
    mcKay
    rename_vertices(::TriangulatedPolygon, ::Vector{<:Integer})
    rename_vertices!(::TriangulatedPolygon, ::Vector{<:Integer})    
    relative_degree(::TriangulatedPolygon, ::Integer, ::Vector{<:Integer})
    relative_degrees(::TriangulatedPolygon, ::Vector{<:Integer}, ::Vector{<:Integer})
```



[^1]: Hartke, S.G., & Radcliffe, A.J. (2008). McKay ’ s Canonical Graph Labeling Algorithm.

# Flip Graphs of Convex Polygons

Flip Graphs, are obtained, by considering triangulations on a fixed set of points as vertices, and connect two vertices, 
if it is possible to get from one triangulation to the other by flipping a single edge.

```@docs
    FlipGraphPlanar
```
`FlipGraphPlanar` implements the AbstractGraph interface from Graphs.jl. It is therefore possible to use it with other packages that work with Graphs.jl. This is very helpfull for plotting the graph.

## Constructors

```@docs
    flipgraph(::TriangulatedPolygon)
    flipgraph_planar
```

## Graph methods

```@docs
    nv(::FlipGraphPlanar)
    ne(::FlipGraphPlanar)
    vertices(::FlipGraphPlanar)
    edges(::FlipGraphPlanar)
    has_vertex(::FlipGraphPlanar,v)
    has_edge(::FlipGraphPlanar,s,d)
    has_edge(::FlipGraphPlanar,::Edge)
    neighbors(::FlipGraphPlanar, ::Integer)
    diameter(::FlipGraphPlanar)
    ```
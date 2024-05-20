
# Flip Graphs of Convex Polygons

Flip Graphs, are obtained, by considering triangulations on a fixed set of points as vertices, and connect two vertices, 
if it is possible to get from one triangulation to the other by flipping a single edge.

```@docs
    FlipGraphPlanar
```
`FlipGraphPlanar` implements the AbstractGraph interface from Graphs.jl. It is therefore possible to use it with other packages that work with Graphs.jl. This is very helpfull for plotting the graph.

```@docs
    flipgraph(::TriangulatedPolygon)
    flipgraph_planar
    nv(::FlipGraphPlanar)
    ne(::FlipGraphPlanar)
    vertices(::FlipGraphPlanar)
    edges(::FlipGraphPlanar)
    has_vertex(::FlipGraphPlanar)
    has_edge(::FlipGraphPlanar)
    neighbors(::FlipGraphPlanar)
    diameter(::FlipGraphPlanar)
```
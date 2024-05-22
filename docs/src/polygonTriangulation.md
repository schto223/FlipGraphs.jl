# Triangulations of Convex Polygons

In order to better understand triangulations, flips and flipgraphs, it can be helpfull to start simple.
If you take any number of points, and iteratively connect them with staight edges that do not cross each other, 
until you can no longer add an edge that does not cross any other edge, what you'll get is a triangulation.

To get from one triangulation of these points to another, you may choose an inner edge, and flip it. 
If you look at any inner edge, the two triangles adjacent to it form a quadrilateral with the edge as one of its diagonals.
By flipping an edge, all we have to do is replace it by the other diagonal. As we don't want to have straight edges that cross each other,
a flip can only be done, if the quadrilateral is convex, and no three of its corners lie on the same line.
As we are interested in taking this theory to closed surfaces, where we will no longer restrict the edges to being straight, we will only consider triangulations of points in convex general position. In this case, a triangulation in the geometric sense is equivalent to a triangulation in the graph sense.

In graph theory, we do not care, if edges are straight or curved. In this case, a triangulation is simply a maximal(in regards to the number of edges) planar graph with a fixed number of points.
Here, we do not care where the points are located, however, in order keep in line with the geometric sense and have a simple visualization,we will consider these points to be the vertices of a convex polygon (i.e. points in convex position). 

```@docs
    TriangulatedPolygon
```
`TriangulatedPolygon` implements the AbstractGraph interface from Graphs.jl. It is therefore possible to use it with other packages that work with Graphs.jl. This is very helpfull for plotting the graph.\\

Vertices are not explicity stored in TriangulatedPolygon. Only the total number of vertices is stored. They are implicitely labeled by the integers from 1 up to the total number of vertices.

Edges are stored as an adjacency list of which vertices are connected to another.

## Constructors 
```@docs
    triangulated_polygon
```

## Graph Methods 

```@docs
    nv(::TriangulatedPolygon)
    ne(::TriangulatedPolygon)
    vertices(::TriangulatedPolygon)
    edges(::TriangulatedPolygon)
    has_vertex(::TriangulatedPolygon, v)
    has_edge(::TriangulatedPolygon, e::Edge)
    has_edge(::TriangulatedPolygon, s, d)
    neighbors(::TriangulatedPolygon, ::Integer)
    degrees(::TriangulatedPolygon)
    is_isomorph(::TriangulatedPolygon, ::TriangulatedPolygon, ::Vector{Vector{Integer}})
    rename_vertices(::TriangulatedPolygon, ::Vector{Integer})
    adjacency_matrix(::TriangulatedPolygon)
```

## Flipping

```@docs    
    is_flippable(::TriangulatedPolygon, ::Integer, ::Integer)
    is_flippable(::TriangulatedPolygon, ::Edge)
    flip!(::TriangulatedPolygon, ::Integer, ::Integer)
    flip!(::TriangulatedPolygon, ::Edge)
    flip(::TriangulatedPolygon, ::Edge)
```
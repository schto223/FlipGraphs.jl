# Triangulations of Convex Polygons

In order to better understand *triangulations*, *flips* and *flip graphs*, it can be helpful to start simple.
If you take any number of points and iteratively connect them with straight edges that do not cross each other 
until you can no longer add an edge that does not cross any other edge, what you'll get is a (geometric) **triangulation**.

To get from one triangulation of a set of points to another, you may choose an inner edge and flip it. 
If you look at any inner edge, the two triangles adjacent to it form a quadrilateral, with the edge as one of its diagonals.
To flip an edge, all we have to do is replace it with the other diagonal. 
As we don't want to have straight edges that cross each other, a flip can only be done if the quadrilateral is convex and no three of its corners lie on the same line.
As we are interested in taking this theory to closed surfaces, where we will no longer have the restriction of edges being straight, we will only consider triangulations of points in convex general position. In this case, a triangulation in the geometric sense is equivalent to a *triangulation of points on the border of a disc* or a *maximal outerplanar graph*.

We do not care where exactly the points are located; however, in order to keep in line with the geometric sense and have a simple visualization, we will consider these points to be the vertices of a convex polygon (i.e. points in convex position). 

## Structures

```@docs
    TriangulatedPolygon
```
`TriangulatedPolygon` implements the `AbstractGraph` interface from Graphs.jl. It is therefore possible to use it with other packages that work with Graphs.jl.

Vertices are not explicitly stored in `TriangulatedPolygon`. 
Only the total number of vertices is stored. They are implicitly labeled by the integers from 1 up to the total number of vertices.\
Edges are stored as an adjacency list.

## Constructors 
```@docs
    triangulated_polygon
```

As an example, the output of `triangulated_polygon(9)` would be a graph that corresponds to the following triangulation of a 9-gon:

```@raw html
<!---![Triangulated 9-gon](assets/triPoly-9.png)-->
<p align="center">
  <img src="assets/triPoly-9.png" width="360px" hspace="20">
</p>
```

## Graph Methods 

The following methods overload some of the main functions from the [Graphs.jl](https://juliagraphs.org/Graphs.jl/stable/) package.

```@docs
    nv(::TriangulatedPolygon)
    ne(::TriangulatedPolygon)
    vertices(::TriangulatedPolygon)
    edges(::TriangulatedPolygon)
    has_vertex(::TriangulatedPolygon, v)
    
    has_edge(::TriangulatedPolygon, e::Edge)
    has_edge(::TriangulatedPolygon, s, d)

    neighbors(::TriangulatedPolygon, ::Integer)
```

If you want to extrude some more information from a `TriangulatedPolygon` object, the following functions might be useful: 

```@docs
    degrees(::TriangulatedPolygon)
    adjacency_matrix(::TriangulatedPolygon)
    diameter(::TriangulatedPolygon)
```

## Flip an edge

```@docs    
    is_flippable(::TriangulatedPolygon, ::Integer, ::Integer)
    is_flippable(::TriangulatedPolygon, ::Edge)
    flip!(::TriangulatedPolygon, ::Integer, ::Integer)    
    flip(::TriangulatedPolygon, ::Integer, ::Integer)
    flip!(::TriangulatedPolygon, ::Edge)
    flip(::TriangulatedPolygon, ::Edge)
    flip_get_edge!
```
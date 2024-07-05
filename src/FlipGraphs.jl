################################################################################
# FlipGraphs.jl
#
#	A Graph-flip package for Julia.
#
# Copyright (C) 2024 Tom Schmit
################################################################################

"""
FlipGraphs

A package for **triangulations** of points and their respective **flip graphs**.

So far, this package includes:\\
- triangulations of **convex polygons** / triangulations of points situated on the border of a disc and their *flip graphs* (with labeled and unlabeled points).\\
- triangulations of points on a **closed surface** and their respective *modular flip graphs* (with labeled and unlabeled points).
"""
module FlipGraphs

import Graphs: 
    Edge, AbstractGraph, SimpleEdge, #diameter,
    has_edge, neighbors, outneighbors, inneighbors, edges, vertices, edgetype, has_vertex, ne, nv, is_directed, src, dst


#general graph functions
export edgetype, has_edge, has_vertex, neighbors, inneighbors, outneighbors, ne, nv
export vertices, edges, is_directed, add_edge!, remove_edge!#, add_vertex!

#graphFunctions
export diameter, adjacency_matrix, invert_permutation, degree

#polygonTriangulations
export TriangulatedPolygon, triangulated_polygon
export degrees
export flip, flip!, flip_get_edge!, is_flippable, edges_flippable, edges_inner, edges_outer, is_inner, is_outer
export mirror!, rotate!, is_identical

#flipGraph_planar
export FlipGraphPlanar, flipgraph, flipgraph_planar
export is_isomorphic, rename_vertices, mcKay

#deltaComplex
export DeltaComplex, deltacomplex, deltacomplex_non_orientable, DualEdge, TriFace
export np, has_point
export get_edge, get_point, get_vertex, get_vertex_id, vertices_id, vertices, edges, points, edges_id, get_edge_id, id
export sides, get_side
export flip!, is_flippable, is_orientable, random_flips!, randomize!, point_degrees, relative_point_degrees
export euler_characteristic, genus, demigenus, diameter_triangulation, diameter_deltaComplex
export adjacency_matrix_triangulation, multi_adjacency_matrix_triangulation, adjacency_matrix_deltacomplex, multi_adjacency_matrix_deltacomplex
export subdivide!, twist_edges!, is_twisted
export rename_edges!, rename_points!, rename_vertices!, is_similar, other_endpoint, quadrilateral_edges

#flipGraph
export FlipGraph, FGVertex, FGVertexCandidate, flipgraph_modular
export mcKay_points, mcKay_vertices, mcKay_edges
export is_isomorphic, vertices_deltacomplex

#general
export matrix_equal, relative_degrees, relative_degree, distances

#exporting
export export_gml


include("generalUtilities.jl")
include("polygonTriangulations.jl")
include("flipGraphPlanar.jl")
include("deltaComplex.jl")
include("flipGraph.jl")
include("exporting.jl")

end # module
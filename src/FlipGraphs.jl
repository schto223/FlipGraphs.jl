################################################################################
# FlipGraphs.jl
#
#	A Graph-flip package for Julia.
#
# Copyright (C) 2024 Tom Schmit
################################################################################

"""
FlipGraphs

A package for triangulations of closed surfaces and their respective FlipGraphs.
"""
module FlipGraphs

import Graphs: 
    Edge, AbstractGraph, SimpleEdge, diameter,
    has_edge, neighbors, outneighbors, inneighbors, edges, vertices, edgetype, has_vertex, ne, nv, is_directed, src, dst


#general graph functions
export edgetype, has_edge, has_vertex, neighbors, inneighbors, outneighbors, ne, nv
export vertices, edges, is_directed, add_edge!, add_vertex!, remove_edge!

#graphFunctions
export diameter, adjacency_matrix, invert_permutation

#polygonTriangulations
export TriangulatedPolygon, triangulated_polygon
export degrees
export flip, flip!, is_flippable

#flipGraph_planar
export FlipGraphPlanar, flipgraph, flipgraph_planar
export is_isomorph, rename_vertices

#deltaComplex
export DeltaComplex, delta_complex, DualEdge, TriFace
export np
export get_edge, get_point, get_vertex, get_vertex_id, vertices_id, vertices, edges, points, edges_id, get_edge_id, id
export flip!, is_flippable, is_orientable, random_flips!, randomize!, point_degrees, relative_point_degrees
export euler_characteristic, genus, demigenus, diameter_triangulation, diameter_deltaComplex
export adjacency_matrix_triangulation, multi_adjacency_matrix_triangulation, adjacency_matrix_deltaComplex
export subdivide!, twist_edges!
export rename_edges!, rename_points!, rename_vertices!, is_similar, other_endpoint

#holeyDeltaComplex
export HoleyDeltaComplex, Hole, Crossing, holey_delta_complex
export num_crossings, edge_crossings, remove_holeloops!, get_crossing

#flipGraph
export FlipGraph, flip_graph
export mcKay_points, mcKay_vertices, mcKay_edges
export is_isomorph, is_isomorph_to


include("generalUtilities.jl")
include("polygonTriangulations.jl")
include("flipGraph_planar.jl")
include("deltaComplex.jl")
include("holeyDeltaComplex.jl")
include("flipGraph.jl")


end # module
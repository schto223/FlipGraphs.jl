################################################################################
# FlipGraphs.jl
#
#	A Graph-flip package for Julia.
#
# Copyright (C) 2024 Tom Schmit
################################################################################

module FlipGraphs

import Graphs: 
    Edge, AbstractGraph, SimpleEdge,  
    has_edge, outneighbors, edges, edgetype, has_vertex, ne, nv, is_directed, vertices, inneighbors, src, dst

import GraphPlot: gplot, cycle_graph
import Cairo, Fontconfig
import Compose: px, PNG, draw

export edgetype, has_edge, has_edge, has_vertex, inneighbors, ne, nv, outneighbors, vertices, is_directed, add_edge!, add_vertex!, remove_edge!



include("planarTriangulations.jl")
include("flipGraph.jl")


end # module
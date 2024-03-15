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
import Compose: px, PNG, draw
import Cairo, Fontconfig

include("planarTriangulations.jl")
include("flipGraph.jl")


end # module
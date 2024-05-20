
using FlipGraphs, Random


HD = holeyDeltaComplex(5, 6)

@time ps = mcKay_triFaces(HD)
@time rename_vertices!(HD, ps[1])

#construct_FlipGraph(HD, 1 ,true)
#@time construct_FlipGraph(HD, 3 ,false)
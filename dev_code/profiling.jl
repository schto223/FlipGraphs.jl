import ProfileView
using FlipGraphs

HD = holeyDeltaComplex(10, 10)
HD2 = deepcopy(HD)
p = mcKay_triFaces(HD; only_one=true)[1]
ProfileView.@profview rename_vertices!(HD, p)
ProfileView.@profview rename_vertices!(HD2, p)

#ProfileView.@ProfileView rename_vertices!(HD2, p)


using ProfileView
#using FlipGraphs

HD = holeyDeltaComplex(5, 6)
HD2 = deepcopy(HD)
p = mcKay_triFaces(HD; only_one=true)
@proofview rename_vertices!(HD, p)
@proofview rename_vertices!(HD2, p)

#ProfileView.@ProfileView rename_vertices!(HD2, p)


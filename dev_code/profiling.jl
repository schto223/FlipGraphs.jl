import ProfileView
using FlipGraphs

D = deltacomplex(1)


ProfileView.@profview flipgraph_modular(D,10)
ProfileView.@profview flipgraph_modular(D,10)


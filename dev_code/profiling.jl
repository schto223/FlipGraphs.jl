import ProfileView
using FlipGraphs

D = deltacomplex(3)


ProfileView.@profview flipgraph_modular(D,2)
ProfileView.@profview flipgraph_modular(D,2)


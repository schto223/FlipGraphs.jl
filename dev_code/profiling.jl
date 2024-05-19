using FlipGraphs
using ProfileView

D = createDeltaComplex(1000,1000)
ProfileView.@profview random_flips!(D,10)
ProfileView.@profview random_flips!(D,10000000)

ProfileView.@profview diameter_triangulation(D)
ProfileView.@profview diameter_triangulation(D)



using FlipGraphs, Random

HD = holey_delta_complex(1,2)
G1 = flip_graph(HD, 0, false)
diameter(G1)
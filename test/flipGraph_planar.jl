@testset "flipGraph_planar" begin

    for n in 3:13
        g = triangulatedPolygon(n)
        G = construct_FlipGraph(g, true)
        # Check if number of triangulations is correct (comparing to https://oeis.org/A000207)
        @test nv(G)== [1, 1, 1, 3, 4, 12, 27, 82, 228, 733, 2282, 7528, 24834, 83898][n-2] 
    end

end
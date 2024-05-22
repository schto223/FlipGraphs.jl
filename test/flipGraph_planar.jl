@testset "flipGraph_planar" begin

    for n in 3:10
        G = flipgraph_planar(n, true)
        # Check if number of triangulations is correct (comparing to https://oeis.org/A000207)
        @test nv(G)== [1, 1, 1, 3, 4, 12, 27, 82, 228, 733, 2282, 7528, 24834, 83898][n-2] 
        @test diameter(G) == [0,0,0,2,2,4,5,6,7,8][n-2]
        
    end

    for n in 3:8
        G = flipgraph_planar(n, false)
        #Comparing to the Catalan numbers: https://oeis.org/A000108
        @test nv(G)== [1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796, 58786, 208012, 742900, 2674440][n-2] 
        @test diameter(G) == [0,1,2,4,5,7,9,11][n-2]
    end
end
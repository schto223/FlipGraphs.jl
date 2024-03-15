@testset "planarTriangulations" begin
	
    for i = 1:10
        g = triangulatedPolygon(i)

	    @test ne(g) == i+(i-3)
    end
end
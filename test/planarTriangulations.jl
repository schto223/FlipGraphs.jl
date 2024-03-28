@testset "planarTriangulations" begin
	
    for i = 1:10
        g = triangulatedPolygon(i)

	    @test ne(g) == i+(i-3)
        @test nv(g) == i

    end

    g = triangulatedPolygon(10)
    for e in edges(g)
        remove_edge!(g, e)  
    end  
    @test ne(g) == 0
    @test nv(g) == 10

    g = triangulatedPolygon(5)
    remove_edge!(g,1,2)
    @test ne(g) == 6
    @test isBiDirectional(g)
    add_edge!(g,1,2)
    @test ne(g) == 7
    @test isBiDirectional(g)

    g = triangulatedPolygon(10)
    FE = filter(e->flippable(g,e) e in edges(g))
    @test length(FE) == 7
    flip!(g, FE[3])
    @test ne(g) = 17
    

end

function isBiDirectional(g::TriGraph)
    for i = 1:ne(g), j in g.adjList[i]
        if !(i in g.adjList[j])
            return false
        end
    end
    return true
end
@testset "flipGraph_planar" begin

    for n in 3:12
        G = flipgraph_planar(n, modular=false)
        #Comparing to the Catalan numbers: https://oeis.org/A000108
        catalan = [1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796, 58786, 208012, 742900, 2674440]
        @test nv(G)== catalan[n-2] 
        @test ne(G)== ((n-3)*catalan[n-2])÷2 
        if n<8
            @test diameter(G) == [0,1,2,4,5,7,9,11][n-2]
        end
    end

    for n in 3:14
        G = flipgraph_planar(n, modular=true)
        # Check if number of triangulations is correct (comparing to https://oeis.org/A000207)
        @test nv(G)== [1, 1, 1, 3, 4, 12, 27, 82, 228, 733, 2282, 7528, 24834, 83898][n-2] 
        @test ne(G)== [0, 0, 0, 2, 4, 20, 66, 248, 868, 3139, 11256, 40832, 148413,543152][n-2]
        if n<9
            @test diameter(G) == [0,0,0,2,2,4,5,6,7,8][n-2]  
        end
    end

    G = flipgraph_planar(8, modular=true)
    
    for e in edges(G)
        @test has_edge(G,e)
    end

    for i in eachindex(G.V)
        for j in eachindex(G.V)
            if i!=j
                @test is_isomorphic(get_vertex(G, i) , get_vertex(G, j)) == false
            end
        end
    end

    @test edgetype(G) == SimpleEdge{Int32}
    @test is_directed(G) == false
    @test is_directed(FlipGraphPlanar) == false

    @test neighbors(G, 5) == inneighbors(G, 5) == outneighbors(G, 5)

    for v in neighbors(G,5)
        @test has_edge(G,5,v)
        @test has_edge(G,v,5)
    end

    for v in vertices(G)
        @test has_vertex(G,v) == true
    end

    @test has_vertex(G, get_vertex(G,3)) == true
    @test has_vertex(G, 3) == true

    g = triangulated_polygon(50)
    p = mcKay(g, only_one=true)[1]
    g2 = rename_vertices(g,p)
    @test is_isomorphic(g,g2)

    g = triangulated_polygon(47)
    perms = mcKay(g)
    g2 = rename_vertices(g,perms[1])
    @test is_isomorphic(g2,g, perms)
    flip!(g2, edges_inner(g2)[4])
    @test is_isomorphic(g2,g, perms) == false
end
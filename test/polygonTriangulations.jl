import Graphs: SimpleEdge, Edge

function isBiDirectional(g::TriangulatedPolygon)
    for i = 1:nv(g), j in g.adjList[i]
        if !(i in g.adjList[j])
            return false
        end
    end
    return true
end

@testset "polygonTriangulations" begin
    for i = 3:20
        g = triangulated_polygon(i)
	    @test ne(g) == i+(i-3)
        @test nv(g) == i
        @test isBiDirectional(g)
    end

    g = triangulated_polygon(5)
    remove_edge!(g,Edge(1,2))
    @test ne(g) == 6
    @test isBiDirectional(g)
    add_edge!(g,1,2)
    @test ne(g) == 7
    @test isBiDirectional(g)

    @testset "flipping" begin
        g = triangulated_polygon(10)
        FE = filter(e -> is_flippable(g,e), edges(g))
        @test length(FE) == 7
        flip!(g, FE[3])
        @test ne(g) == 17
        FE = filter(e->is_flippable(g,e), edges(g))
        flip!(g, FE[3].src, FE[3].dst)

        for e in edges(g)
            @test has_edge(g,e) == true
            if is_flippable(g,e)
                gg = flip(g,e)
                @test nv(gg) == nv(g)
            end
        end
    end

    g = triangulated_polygon(10)
    #vertices
    for i in vertices(g)
        @test has_vertex(g,i) == true
    end
    @test has_vertex(g, nv(g)+1) == false
    @test has_vertex(g, 0) == false

    A = adjacency_matrix(g)
    @test size(A) == (10,10)
    @test ((A.==1) + (A.==0)) == ones(Int, 10, 10)

    @test edgetype(g) == SimpleEdge{eltype(g)}
    @test is_directed(g) == false
    @test is_directed(TriangulatedPolygon) == false

    @test outneighbors(g,4) == inneighbors(g,4) == neighbors(g,4)

    #edges_outer, edges_inner
    g = triangulated_polygon(3)
    @test edges_outer(g) == edges(g)
    @test isempty(edges_inner(g)) == true
    g = triangulated_polygon(27)
    @test isempty(intersect(edges_outer(g), edges_inner(g)))
    @test ne(g) == length(edges_outer(g)) + length(edges_inner(g))
    g = triangulated_polygon(7)
    for e in edges_inner(g)
        @test is_inner(g, e) == true
    end
    for e in edges_outer(g)
        @test is_outer(g, e) == true
    end

    #is_identical
    g1 = triangulated_polygon(30)
    g2 = triangulated_polygon(31)
    @test is_identical(g1,g2) == false
    g2 = flip(g1, edges_inner(g1)[5])
    @test is_identical(g1,g2) == false

    @testset "mcKay" begin
        g = triangulated_polygon(100)
        g_copy= deepcopy(g)
        p = mcKay(g, only_one=true)[1]
        p_inv = invert_permutation(p)
        rename_vertices!(g, p)
        @test g.adjList == rename_vertices(g_copy, p).adjList
        rename_vertices!(g, p_inv)
        @test g.adjList == g_copy.adjList
    end
    
    @testset "symmetry" begin
        for n in 3:10
            g = triangulated_polygon(n)
            g2 = deepcopy(g)
            for i in 1:n
                rotate!(g,1)
                @test is_isomorphic(g, g2)
            end
            @test is_identical(g, g2) == true
            @test is_identical(rotate!(g,n), g2) == true
            mirror!(g)
            @test is_isomorphic(g, g2)
        end
    end

end
using Random

@testset "deltaComplex" begin

    @testset "Sphere" begin
        for p = 3:8
            D = createDeltaComplex(0, p)
            @test nv(D) == 2*(p-2)  
            @test ne(D) == 3*(p-2) 
            @test np(D) == p
            @test is_orientable(D) == true
            @test euler_characteristic(D) == 2
            @test genus(D) == 0
        end
    end

    @testset "Orientable surfaces" begin
        for g = 1:7
            for p = 1:7
                D = createDeltaComplex(g, p)
                @test genus(D) == g
                @test np(D) == p
                @test is_orientable(D) == true
                @test euler_characteristic(D) == 2-2*g         
            end
        end

        #test flip
        D = createDeltaComplex(3,7)
        @test sum(point_degrees(D)) == 2*ne(D)
        e = get_edge(D,1)
        e_copy = deepcopy(e)
        t1,t2 = vertices(e)
        T1 = deepcopy(get_vertex(D, t1))
        T2 = deepcopy(get_vertex(D, t2))
        flip!(D,e)
        flip!(D,e)
        @test genus(D) == 3
        @test np(D) == 7
        @test is_orientable(D) == true
        @test euler_characteristic(D) == -4  
        flip!(D,e)
        flip!(D,e)        
        @test all(e_copy.triangles.==e.triangles) && all(e_copy.sides==e.sides) && e.is_twisted==e_copy.is_twisted

        #test createDeltaComplex 
        D = createDeltaComplex([1,2,3,-3,-1,-2]) # torus with an added point that has only one outgoing edge
        @test sum(point_degrees(D)) == 2*ne(D)
        @test np(D)==2
        @test genus(D)==1
        e = filter(e-> 2 ∈ points(D,e) , edges(D))[1]
        @test is_flippable(e) == false
    end

    @testset "non-orientable surfaces" begin
        @testset "kleinBottle" begin
            D = createDeltaComplex([1,2,-1,2]) #klein Bottle
            @test is_orientable(D) == false
            @test np(D) == 1 && nv(D) == 2 && ne(D) == 3
            @test demigenus(D) == 2
            e = filter(e-> !e.is_twisted , edges(D))[1]
            flip!(D,e) 
            @test is_orientable(D) == false
            @test np(D) == 1 && nv(D) == 2 && ne(D) == 3
            @test demigenus(D) == 2
            e = filter(e-> e.is_twisted , edges(D))[1]
            flip!(D,e) 
            @test is_orientable(D) == false
            @test np(D) == 1 && nv(D) == 2 && ne(D) == 3
            @test demigenus(D) == 2
        end

        @testset "projective plane" begin
            D = createDeltaComplex([1,2,1,2]) # projective plane
            @test is_orientable(D) == false
            @test np(D) == 2 && nv(D) == 2 && ne(D) == 3
            @test demigenus(D) == 1
            e = filter(e-> !e.is_twisted , edges(D))[1]
            flip!(D,e) 
            @test is_orientable(D) == false
            @test np(D) == 2 && nv(D) == 2 && ne(D) == 3
            @test demigenus(D) == 1
            e = filter(e-> e.is_twisted , edges(D))[1]
            flip!(D,e) 
            @test is_orientable(D) == false
            @test np(D) == 2 && nv(D) == 2 && ne(D) == 3
            @test demigenus(D) == 1
        end

        @testset "projective plane * kleinBottle" begin
            D = createDeltaComplex([1,2,1,2,3,4,-3,4]) # projective plane glued to klein Bottle
            @test is_orientable(D) == false
            @test np(D) == 2 && nv(D) == 6 && ne(D) == 9
            @test sum(point_degrees(D)) == 2*ne(D)
            @test demigenus(D) == 3
            e = filter(e-> !e.is_twisted , edges(D))[1]
            flip!(D,e) 
            @test is_orientable(D) == false
            @test np(D) == 2 && nv(D) == 6 && ne(D) == 9
            @test demigenus(D) == 3
            e = filter(e-> e.is_twisted , edges(D))[1]
            flip!(D,e) 
            @test is_orientable(D) == false
            @test np(D) == 2 && nv(D) == 6 && ne(D) == 9
            @test demigenus(D) == 3
        end
    end

    @testset "Errors" begin
        @test_throws ArgumentError createDeltaComplex([1])
        @test_throws ArgumentError createDeltaComplex([1,2,1,3])
        @test_throws ArgumentError createDeltaComplex([1,2,1,-2,1])

        D = createDeltaComplex(4)
        @test_throws ArgumentError demigenus(D)
        @test_throws ArgumentError createDeltaComplex(0)
        D = createDeltaComplex(0,3)
        @test_throws ArgumentError demigenus(D)

        D = createDeltaComplex([1,2,1,2])
        @test_throws ArgumentError genus(D)
        D = createDeltaComplex([1,2,-1,2])
        @test_throws ArgumentError genus(D)
    end

    @testset "diameter" begin
        D = createDeltaComplex(10,20)
        @test 1 <= diameter_triangulation(D) <= np(D)-1
        @test 1 <= diameter_deltaComplex(D) <= nv(D)-1
        Random.seed!(1234)
        a = rand(1:ne(D), 1000)
        for i in eachindex(a)
            if is_flippable(D,a[i])
                flip!(D,a[i])
            end
        end 
        @test 1 <= diameter_triangulation(D) <= np(D)-1
        @test 1 <= diameter_deltaComplex(D) <= nv(D)-1
    end

    @testset "Random flipping" begin
        D1 = createDeltaComplex(10,50)
        D2 = createDeltaComplex(10,50)

        randomize!(D1, 100000,10000,10)
        @test nv(D1) == nv(D2)
        @test ne(D1) == ne(D2)
        @test np(D1) == np(D2)
        @test sum(point_degrees(D1)) == 2*ne(D1)
    end

end


    

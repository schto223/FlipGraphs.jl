@testset "deltaComplex" begin

    @testset "Sphere" begin
        for p=3:10
            D = createDeltaComplex(0,p)
            @test nv(D) == 2*(p-2)  
            @test ne(D) == 3*(p-2) 
            @test np(D) == p
            @test is_orientable(D) == true
            @test euler_characteristic(D) == 2
            @test genus(D) == 0
        end
    end

    @testset "Orientable surfaces" begin
        for g = 1:8
            for p = 1:8
                D = createDeltaComplex(g,p)
                @test genus(D) == g
                @test np(D) == p
                @test is_orientable(D) == true
                @test euler_characteristic(D) == 2-2*g         
            end
        end

        D = createDeltaComplex(3,7)
        e = get_edge(D,1)
        e_copy = deepcopy(e)
        t1,t2 = trifaces(e)
        T1 = deepcopy(get_vertex(D, t1))
        T2 = deepcopy(get_vertex(D, t2))
        flip!(D,e)
        flip!(D,e)
        #@test e_copy == e
        @test genus(D) == 3
        @test np(D) == 7
        @test is_orientable(D) == true
        @test euler_characteristic(D) == -4  

        flip!(D,e)
        flip!(D,e)
        #@test T1 == get_vertex(D, t1)
        #@test T2 == get_vertex(D, t2)
        
    end

end


    

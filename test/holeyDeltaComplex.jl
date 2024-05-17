using Random

function holes_intact(HD::HoleyDeltaComplex)
    for H in HD.holes
        c_first = get_crossing(H)
        c = c_first
        skipfirst = true
        while c != c_first || skipfirst
            skipfirst = false
            if c.previous.next != c || c.next.previous != c
                return false
            end
            d_prev = get_edge(HD, c.previous.edge_id)
            d = get_edge(HD, c.edge_id)
            if !(d.triangles[1] in vertices(d_prev) || d.triangles[2] in vertices(d_prev))
                return false
            end
            c = c.next
        end
    end
    return true
end  

@testset "HoleyDeltaComplex" begin

    #@testset "Sphere" begin
    #    for p=3:10
    #        D = createDeltaComplex(0, p)
    #        @test nv(D) == 2*(p-2)  
    #        @test ne(D) == 3*(p-2) 
    #        @test np(D) == p
    #        @test euler_characteristic(D) == 2
    #        @test genus(D) == 0
    #    end
    #end

    @testset "Orientable surfaces" begin
        for g = 1:8
            for p = 1:8
                HD = holeyDeltaComplex(g, p)
                @test genus(HD) == g
                @test np(HD) == p
                @test euler_characteristic(HD) == 2-2*g    
                @test sum(H->num_crossings(H), HD.holes) == length(reduce(vcat, HD.edge_crossings)) 
                @test holes_intact(HD)   
            end
        end

        #test flip
        HD = holeyDeltaComplex(3,7)
        @test sum(point_degrees(HD)) == 2*ne(HD)
        e = get_edge(HD,1)
        e_copy = deepcopy(e)
        t1,t2 = vertices(e)
        T1 = deepcopy(get_vertex(HD, t1))
        T2 = deepcopy(get_vertex(HD, t2))
        flip!(HD,e)
        @test genus(HD) == 3
        @test np(HD) == 7
        @test euler_characteristic(HD) == -4  
        @test sum(H->num_crossings(H), HD.holes) == length(reduce(vcat, HD.edge_crossings))  
        @test holes_intact(HD)   
        flip!(HD,e, false)        
        @test all(e_copy.triangles.==e.triangles) && all(e_copy.sides==e.sides) && e.is_twisted==e_copy.is_twisted
        @test sum(H->num_crossings(H), HD.holes) == length(reduce(vcat, HD.edge_crossings))  
        @test holes_intact(HD)   
    end 

    @testset "diameter & reversability" begin 
        HD = holeyDeltaComplex(10,20)
        @test 1 <= diameter_triangulation(HD) <= np(HD)-1
        @test 1 <= diameter_deltaComplex(HD) <= nv(HD)-1
        Random.seed!(1234)
        a = rand(1:ne(HD), 1000)
        dir = rand(Bool, 1000)
        for i in eachindex(a)
            flip!(HD, a[i], dir[i])
        end 
        @test 1 <= diameter_triangulation(HD) <= np(HD)-1
        @test 1 <= diameter_deltaComplex(HD) <= nv(HD)-1
        @test sum(H->num_crossings(H), HD.holes) == length(reduce(vcat, HD.edge_crossings))  
        @test holes_intact(HD)   

        for i in ne(HD):-1:1
            flip!(HD, a[i], !dir[i])
        end
    end
end


  

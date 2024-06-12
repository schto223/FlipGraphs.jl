using Random

@testset "generalUtilities" begin
    Random.seed!(1234)
    A = rand([0,0,0,0,0,0,0,0,0,1], 1000, 1000)
    A = transpose(A) + A
    for i in eachindex(A)
        if A[i] ==2
            A[i] = 1
        end
    end
    D = distances(A)
    @test all(D.>=0)
    @test maximum(D)<1000
end
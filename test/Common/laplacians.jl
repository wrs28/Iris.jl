using Iris
using Test

@testset "Laplacians" begin
    @testset "1D" begin
        @testset "Cartesian" begin
            lat = Lattices.Cartesian(.1)
            lap = @test_nowarn Laplacian{Iris.Common.Laplacians.Symmetric}(lat,ones(ComplexF64,10),ones(ComplexF64,11))
            for i ∈ 1:10 @test lap.l0[i,i] ≈ -200 end
            for i ∈ 1:9 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 2:10 @test lap.l0[i,i-1] ≈ 100 end
        end
    end
    @testset "2D" begin
        @testset "Cartesian" begin
            lat = Cartesian(.1,.2)
            lap = @test_nowarn Laplacian{Iris.Common.Laplacians.Symmetric}(lat,ones(ComplexF64,10),ones(ComplexF64,11),ones(ComplexF64,11),ones(ComplexF64,12))
            for i ∈ 1:10*11 @test lap.l0[i,i] ≈ -250 end
            for i ∈ 1:9 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 1:9 @test lap.l0[i,i+10] ≈ 25 end
            for i ∈ 11:19 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 11:19 @test lap.l0[i,i+10] ≈ 25 end
            for i ∈ 21:29 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 21:29 @test lap.l0[i,i+10] ≈ 25 end
            for i ∈ 31:39 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 31:39 @test lap.l0[i,i+10] ≈ 25 end
            for i ∈ 41:49 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 41:49 @test lap.l0[i,i+10] ≈ 25 end
            for i ∈ 51:59 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 51:59 @test lap.l0[i,i+10] ≈ 25 end
            for i ∈ 61:69 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 61:69 @test lap.l0[i,i+10] ≈ 25 end
            for i ∈ 71:79 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 71:79 @test lap.l0[i,i+10] ≈ 25 end
            for i ∈ 81:89 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 81:89 @test lap.l0[i,i+10] ≈ 25 end
            for i ∈ 91:99 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 91:99 @test lap.l0[i,i+10] ≈ 25 end
            for i ∈ 101:109 @test lap.l0[i,i+1] ≈ 100 end
            for i ∈ 2:10 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 12:20 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 12:20 @test lap.l0[i,i-10] ≈ 25 end
            for i ∈ 22:30 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 22:30 @test lap.l0[i,i-10] ≈ 25 end
            for i ∈ 32:40 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 32:40 @test lap.l0[i,i-10] ≈ 25 end
            for i ∈ 42:50 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 42:50 @test lap.l0[i,i-10] ≈ 25 end
            for i ∈ 52:60 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 52:60 @test lap.l0[i,i-10] ≈ 25 end
            for i ∈ 62:70 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 62:70 @test lap.l0[i,i-10] ≈ 25 end
            for i ∈ 72:80 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 72:80 @test lap.l0[i,i-10] ≈ 25 end
            for i ∈ 82:90 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 82:90 @test lap.l0[i,i-10] ≈ 25 end
            for i ∈ 92:100 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 92:100 @test lap.l0[i,i-10] ≈ 25 end
            for i ∈ 102:110 @test lap.l0[i,i-1] ≈ 100 end
            for i ∈ 102:110 @test lap.l0[i,i-10] ≈ 25 end
            lap.l0[1,1] = 1im
            @test_nowarn conj!(lap)
            @test lap.l0[1,1] ≈ -1im
        end
    end
end

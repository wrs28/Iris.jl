using Iris
using Test
using LinearAlgebra
using SparseArrays

@testset "LU Factorizations" begin
    A = I + sprand(1000,1000,.2)
    y = rand(1000)
    alu = @test_nowarn lu(A,USolver)
    @test norm(A*(alu\y)-y) ≤ 1e-12
    # alu = @test_skip lu(A,PSolver)
    # @test_skip norm(A*(alu\y)-y) ≤ 1e-12
    # alu = @test_skip lu(A,MSolver)
    # @test_skip norm(A*(alu\y)-y) ≤ 1e-12
end

using Iris
using Test

@testset "PumpFunctions" begin
        F = @test_nowarn Iris.Common.PumpFunctions.PiecewiseConstant()
        @test F(Point(0,0)) ≈ 1
        F = @test_nowarn Iris.Common.PumpFunctions.PiecewiseConstant(.1)
        @test F(Point(0,0)) ≈ .1
end

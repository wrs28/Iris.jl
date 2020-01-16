using Iris
using Test

@testset "DielectricFunctions" begin
        ε = @test_nowarn Iris.Common.DielectricFunctions.PiecewiseConstant(1)
        @test ε(Point(0,0)) ≈ complex(1,0)^2
        ε = @test_nowarn Iris.Common.DielectricFunctions.PiecewiseConstant(1,1)
        @test ε(Point(0,0)) ≈ complex(1,1)^2
        ε = @test_nowarn Iris.Common.DielectricFunctions.PiecewiseConstant(1+.2im)
        @test ε(Point(0,0)) ≈ complex(1,.2)^2
        @test_nowarn ε.n1 = .1
        @test ε(Point(0,0)) ≈ complex(.1,.2)^2
        @test_nowarn ε.n₁ = .2
        @test ε(Point(0,0)) ≈ complex(.2,.2)^2
        @test_nowarn ε.n2 = .3
        @test ε(Point(0,0)) ≈ complex(.2,.3)^2
        @test_nowarn ε.n₂ = .4
        @test ε(Point(0,1000)) ≈ complex(.2,.4)^2
        @test_nowarn ε.ε = complex(2,4)^2
        @test ε(Point(0,1000)) ≈ complex(2,4)^2
        @test ε.n1 ≈ 2
        @test ε.n₁ ≈ 2
        @test ε.n2 ≈ 4
        @test ε.n₂ ≈ 4
end

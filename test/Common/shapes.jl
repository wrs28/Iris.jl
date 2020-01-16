using Iris
using Test

@testset "Shapes" begin
    @testset "1D" begin
        @testset "Interval" begin
            i = Interval(0,1; reference=:left)
            @test i.origin.x≈0
            i = Interval(0,1; reference=:right)
            @test i.origin.x≈1
            i = Interval(0,1; reference=:center)
            @test i.origin.x≈.5
            @test i(.5)
            @test !i(-.5)
            @test !i(1.5)
            @test ndims(i)==1
            @test Iris.Common.nsides(i) == 2
        end
    end
    @testset "2D" begin
        @testset "Circle" begin
            c = Circle(1,0,0; reference=:bottom)
            @test c.origin≈Point(0,1)
            c = Circle(1,0,0; reference=:top)
            @test c.origin≈Point(0,-1)
            c = Circle(5,Point(0,0); ϕ=π/3, reference=:left)
            @test c.origin≈Point(5/2,5sqrt(3)/2)
            c = Circle(10,0,0; reference=:right)
            @test c.origin≈Point(-10,0)
            @test !c(Point(.1,0))
            @test c(Point(-.1,0))
            @test c(Point(-10,10-2eps(10.)))
            @test !c(Point(-10,-10-2eps(10.)))
        end
        @testset "Square" begin
            s = Square(1,0,0; reference=:bottom)
            @test s.origin≈Point(0,1/2)
            s = Square(1,0,0; reference=:top)
            @test s.origin≈Point(0,-1/2)
            s = Square(2,Point(0,0); ϕ=π/3, reference=:left)
            @test s.origin≈Point(1/2,1sqrt(3)/2)
            s = Square(10,0,0; reference=:right)
            @test s.origin≈Point(-10/2,0)
            @test !s(Point(.1,0))
            @test s(Point(-.1,0))
            @test s(Point(-5,5-2eps(5.)))
        end
        @testset "Rectangle" begin
            r = Rectangle(1,.5,0,0; reference=:bottom)
            @test r.origin≈Point(0,.25)
            r = Rectangle(1,.5,0,0; reference=:top)
            @test r.origin≈Point(0,-.25)
            r = Rectangle(2,1,Point(0,0); ϕ=π/3, reference=:left)
            @test r.origin≈Point(1/2,sqrt(3)/2)
            r = Rectangle(10,1,0,0; reference=:right)
            @test r.origin≈Point(-10/2,0)
            @test !r(Point(.1,0))
            @test r(Point(-.1,0))
            @test r(Point(-1,.5-2eps(.5)))
            @test !r(Point(-1,-1-2eps(1.)))
        end
    end
end

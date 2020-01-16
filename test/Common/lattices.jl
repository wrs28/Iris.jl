using Iris
using Test

@testset "Lattices" begin
    @testset "1D" begin
        @testset "Cartesian" begin
            lat = @test_nowarn Lattices.Cartesian(.313)
            @test lat[latticeindex(lat,π)] ≈ Point(π)
            @test all(lat[0:2] .≈ [Point(0),Point(.313),Point(.626)])
            lat = @test_nowarn Lattices.Cartesian(.1; origin = Point(2))
            @test lat[latticeindex(lat,π)] ≈ Point(π)
            @test_throws ArgumentError Lattices.Cartesian(.1; origin = Point(2), angles = .1)
        end
        @testset "Polar" begin
            lat = @test_nowarn Lattices.Polar(.313)
            @test lat[latticeindex(lat,π)] ≈ Point(π)
            @test all(lat[0:2] .≈ [Point(lat.dr/2),Point(3lat.dr/2),Point(5lat.dr/2)])
            lat = @test_nowarn Lattices.Polar(.1; origin = Point(2))
            @test lat[latticeindex(lat,π)] ≈ Point(π)
            @test_throws ArgumentError Lattices.Polar(.1; origin = Point(2), angles = .1)
        end
        @testset "Spherical" begin
            lat = @test_nowarn Lattices.Spherical(.313)
            @test lat[latticeindex(lat,π)] ≈ Point(π)
            @test all(lat[0:2] .≈ [Point(lat.dr/2),Point(3lat.dr/2),Point(5lat.dr/2)])
            lat = @test_nowarn Lattices.Spherical(.1; origin = Point(2))
            @test lat[latticeindex(lat,π)] ≈ Point(π)
            @test_throws ArgumentError Lattices.Spherical(.1; origin = Point(2), angles = .1)
        end
    end
    @testset "2D" begin
        @testset "Cartesian" begin
            lat = @test_nowarn Lattices.Cartesian(.313,.213)
            @test lat[latticeindex(lat,Point(π,2))] ≈ Point(π,2)
            @test all(lat[0:2,0:2] .≈ [Point(0,0),Point(.313,.213),Point(.626,.426)])
            @test all([lat[i,j] for i ∈ 0:1, j ∈ 0:1] .≈ [Point(0,0) Point(0,.213);Point(.313,0) Point(.313,.213)])
            lat = @test_nowarn Lattices.Cartesian(.1,.2; origin = Point(2,3))
            @test lat[latticeindex(lat,Point(π,sqrt(5)))] ≈ Point(π,sqrt(5))
            @test_nowarn Lattices.Cartesian(.1,.2; origin = Point(2,3), angles = .1)
        end
        @testset "Polar" begin
            lat = @test_nowarn Lattices.Polar(.1,.1)
            @test lat[latticeindex(lat,Point{Polar}(1,π))] ≈ Point(-1,0)
            @test all(lat[0:2,0:2] .≈ [Point{Polar}(.05,0),Point{Polar}(0.15,.1),Point{Polar}(.25,.2)])
            lat = @test_nowarn Lattices.Polar(.1; origin = Point(2))
            @test lat[latticeindex(lat,π)] ≈ Point(π)
            @test_throws ArgumentError Lattices.Polar(.1; origin = Point(2), angles = .1)
        end
        @testset "Spherical" begin
            lat = @test_nowarn Lattices.Spherical(1,π/4)
            @test lat[latticeindex(lat,1,π)] ≈ Point(0,-1)
            @test all(lat[0:2,0:2] .≈ [Point(0,1/2),Point(1.5/sqrt(2),1.5/sqrt(2)),Point(2.5,0)])
            lat = @test_nowarn Lattices.Spherical(.1,.2; origin = Point(2,1))
            @test lat[latticeindex(lat,π,π)] ≈ Point(0,-π)
            @test_throws ArgumentError Lattices.Spherical(.1; origin = Point(2), angles = .1)
        end
    end
end

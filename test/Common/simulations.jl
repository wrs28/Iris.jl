using Iris
using Test

@testset "Simulations" begin
    @testset "1D" begin
        @testset "Cartesian" begin
            cavity = Interval(-.5,.5)
            pump = Interval(-.5,0)
            bnd = Boundary(Interval(-1,1),PML)
            dom = NondispersiveDomain(cavity,2)
            fdom = DispersiveDomain(pump,TwoLevelSystem(10,.1,2))
            lat = Cartesian(.05)
            latdom = LatticeDomain(bnd,lat)
            sim = @test_nowarn Simulation(10,dom,latdom,fdom)
            @test_nowarn smooth!(sim)
            dom.n=3
            fdom.F=2
            @test_nowarn update!(sim)
            @test_nowarn smooth!(sim)
            @test_nowarn update_dielectric!(sim,rand(ComplexF64,40))
            @test_nowarn update_pump!(sim,rand(Float64,40))
            @test_nowarn smooth!(sim)
        end
    end
    @testset "2D" begin
        @testset "Cartesian" begin
            cavity = Circle(.5,0,0)
            pump = Square(.3,0,0)
            bnd = Boundary(Rectangle(1,.7,0,0),PML)
            dom = NondispersiveDomain(cavity,2)
            fdom = DispersiveDomain(pump,TwoLevelSystem(10,.1,2))
            lat = Cartesian(.025,.035)
            latdom = LatticeDomain(bnd,lat)
            sim = @test_nowarn Simulation(10,dom,latdom,fdom)
            # @test_nowarn smooth!(sim)
            dom.n=3
            fdom.F=2
            @test_nowarn update!(sim)
            # @test_nowarn smooth!(sim)
            @test_nowarn update_dielectric!(sim,rand(ComplexF64,40))
            @test_nowarn update_pump!(sim,rand(Float64,40))
            # @test_nowarn smooth!(sim)
        end
    end
end

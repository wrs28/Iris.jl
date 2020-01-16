using Iris
using Test

@testset "Helmholtz Operators" begin
    @testset "1D" begin
        cavity = Interval(-.5,.5)
        pump = Interval(-.5,0)
        bnd = Boundary(Interval(-1,1),noBL,DirichletBC)
        dom = NondispersiveDomain(cavity,2)
        fdom = DispersiveDomain(pump,TwoLevelSystem(10,.1,2))
        lat = Cartesian(.05)
        latdom = LatticeDomain(bnd,lat)
        sim = Simulation(10,dom,latdom,fdom)
        h = @test_nowarn Helmholtz(sim; m=2)
        d = @test_nowarn h(0)
        @test d[1,1] ≈ -1200
        @test d[1,2] ≈ 400
        @test d[2,1] ≈ 400
        @test d[2,2] ≈ -800
        bnd = Boundary(Interval(-1,1),noBL,NeumannBC)
        latdom = LatticeDomain(bnd,lat)
        sim = Simulation(10,dom,latdom,fdom)
        h = @test_nowarn Helmholtz(sim; m=2)
        d = @test_nowarn h(0)
        @test d[1,1] ≈ -400
        @test d[1,2] ≈ 400
        @test d[2,1] ≈ 400
        @test d[2,2] ≈ -800
        bnd = Boundary(Interval(-1,1),noBL,MatchedBC{1}(in=[1]),MatchedBC{2}(in=[1]))
        latdom = LatticeDomain(bnd,lat)
        sim = Simulation(10,dom,latdom,fdom)
        h = @test_nowarn Helmholtz(sim; m=2)
        d = @test_nowarn h(20)
    end
    @testset "2D" begin
        cavity = Circle(.4,0,0)
        pump = Square(.3,0,0)
        bnd = Boundary(Rectangle(1,.9,0,0),PML)
        dom = NondispersiveDomain(cavity,2)
        fdom = DispersiveDomain(pump,TwoLevelSystem(10,.1,2))
        lat = Cartesian(.025,.035)
        latdom = LatticeDomain(bnd,lat)
        sim = Simulation(10,dom,latdom,fdom)
        h = @test_nowarn Helmholtz(sim; m=1)
    end
end

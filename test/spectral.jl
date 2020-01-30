using Iris
using Test

@testset "Linear Eigenvalue Problem" begin
    @testset "Helmholtz 1D" begin
        @testset "PML" begin
            cavity = Interval(-.5,.5)
            pump = Interval(-.5,0)
            bnd = Boundary(Interval(-1,1),PML)
            dom = NondispersiveDomain(cavity,3)
            lat = Cartesian(.001)
            latdom = LatticeDomain(bnd,lat)
            sim = Simulation(10,dom,latdom)
            lep = @test_nowarn HelmholtzLEP(sim)
            ω, ψ = @test_nowarn helmholtzeigen(lep,10;nev=4)
            sort!(ω; by=real)
            @test isapprox(ω[2] - ω[1], π/3; atol=1e-3)
            @test isapprox(ω[3] - ω[2], π/3; atol=1e-3)
            @test isapprox(ω[4] - ω[3], π/3; atol=1e-3)
        end
        @testset "cPML" begin
            cavity = Interval(-.5,.5)
            pump = Interval(-.5,0)
            bnd = Boundary(Interval(-1,1),cPML)
            dom = NondispersiveDomain(cavity,3)
            lat = Cartesian(.001)
            latdom = LatticeDomain(bnd,lat)
            sim = Simulation(10,dom,latdom)
            lep = @test_nowarn HelmholtzLEP(sim)
            ω, ψ = @test_nowarn helmholtzeigen(lep,10;nev=4)
            sort!(ω; by=real)
            @test isapprox(ω[2] - ω[1], π/3; atol=1e-3)
            @test isapprox(ω[3] - ω[2], π/3; atol=1e-3)
            @test isapprox(ω[4] - ω[3], π/3; atol=1e-3)
        end
    end
    @testset "Helmholtz 2D" begin
        @testset "Dirichlet" begin
            bnd = Boundary(Square(π,0,0),DirichletBC)
            lat = Cartesian(.004,.004)
            latdom = LatticeDomain(bnd,lat)
            sim = Simulation(1,latdom)
            lep = @test_nowarn HelmholtzLEP(sim)
            ω, ψ = @test_nowarn helmholtzeigen(lep,1;nev=10)
            sort!(ω; by=real)
            @test isapprox(ω[1]^2, 1^2 + 1^2; atol = 3e-2)
            @test isapprox(ω[2]^2, 2^2 + 1^2; atol = 3e-2)
            @test isapprox(ω[3]^2, 1^2 + 2^2; atol = 3e-2)
            @test isapprox(ω[4]^2, 2^2 + 2^2; atol = 3e-2)
            @test isapprox(ω[5]^2, 1^2 + 3^2; atol = 3e-2)
            @test isapprox(ω[6]^2, 3^2 + 1^2; atol = 3e-2)
            @test isapprox(ω[7]^2, 2^2 + 3^2; atol = 3e-2)
            @test isapprox(ω[8]^2, 3^2 + 2^2; atol = 3e-2)
            @test isapprox(ω[9]^2, 1^2 + 4^2; atol = 3e-2)
            @test isapprox(ω[10]^2, 4^2 + 1^2; atol = 3e-2)
        end
    end
end
@testset "CF Eigenvalue Problem" begin
    @testset "Helmholtz 1D" begin
        @testset "PML" begin
            cavity = Interval(-.5,.5)
            pump = Interval(-.5,0)
            bnd = Boundary(Interval(-1,1),PML)
            dom = NondispersiveDomain(cavity,3)
            pump = DispersiveDomain(Interval(-.8,.8),TwoLevelSystem())
            lat = Cartesian(.001)
            latdom = LatticeDomain(bnd,lat)
            sim = Simulation(10,dom,latdom,pump)
            cf = @test_nowarn HelmholtzCF(sim)
            η , u = @test_nowarn helmholtzeigen(cf,10;nev=4)
        end
        @testset "cPML" begin
            cavity = Interval(-.5,.5)
            pump = Interval(-.5,0)
            bnd = Boundary(Interval(-1,1),cPML)
            dom = NondispersiveDomain(cavity,3)
            pump = DispersiveDomain(Interval(-.8,.8),TwoLevelSystem())
            lat = Cartesian(.001)
            latdom = LatticeDomain(bnd,lat)
            sim = Simulation(10,dom,latdom,pump)
            cf = @test_nowarn HelmholtzCF(sim)
            η , u = @test_nowarn helmholtzeigen(cf,10;nev=4)
        end
    end
end
@testset "Nonlinear Eigenvalue Problem" begin
    @testset "Helmholtz 1D" begin
        @testset "PML" begin
            cavity = Interval(-.5,.5)
            pump = Interval(-.5,0)
            bnd = Boundary(Interval(-1,1),PML)
            dom = NondispersiveDomain(cavity,3)
            fdom = DispersiveDomain(pump,TwoLevelSystem(10,.1,2))
            lat = Cartesian(.001)
            latdom = LatticeDomain(bnd,lat)
            sim = Simulation(10,dom,latdom,fdom)
            nep = @test_nowarn HelmholtzNEP(sim)
            @test_nowarn nep(10)
            ω, ψ = @test_nowarn helmholtzeigen(nep,10;neigs=4)
        end
        @testset "cPML" begin
            cavity = Interval(-.5,.5)
            pump = Interval(-.5,0)
            bnd = Boundary(Interval(-1,1),cPML)
            dom = NondispersiveDomain(cavity,3)
            fdom = DispersiveDomain(pump,TwoLevelSystem(10,.1,2))
            lat = Cartesian(.001)
            latdom = LatticeDomain(bnd,lat)
            sim = Simulation(10,dom,latdom,fdom)
            nep = @test_nowarn HelmholtzNEP(sim)
            @test_nowarn nep(10)
            ω, ψ = @test_nowarn helmholtzeigen(nep,10;neigs=4)
        end
        @testset "MathcedBC sans dispersion" begin
            cavity = Interval(-.5,.5)
            pump = Interval(-.5,0)
            bnd = Boundary(Interval(-1,1),MatchedBC{1}(out=1),MatchedBC{2}(out=1))
            dom = NondispersiveDomain(cavity,3)
            lat = Cartesian(.001)
            latdom = LatticeDomain(bnd,lat)
            sim = Simulation(10,dom,latdom)
            nep = @test_nowarn HelmholtzNEP(sim)
            @test_nowarn nep(10)
            ω, ψ = @test_nowarn helmholtzeigen(nep,10;neigs=4,maxit=50)
            sort!(ω; by=real)
            @test isapprox(ω[2] - ω[1], π/3; atol=1e-3)
            @test isapprox(ω[3] - ω[2], π/3; atol=1e-3)
            @test isapprox(ω[4] - ω[3], π/3; atol=1e-3)
        end
        @testset "MathcedBC with complex dispersion" begin
            cavity = Interval(-.5,.5)
            bnd = Boundary(Interval(-1,1),MatchedBC{1}(out=1),MatchedBC{2}(out=1))
            dom = NondispersiveDomain(cavity,3)
            fdom1 = DispersiveDomain(Interval(-.5,-.1),TwoLevelSystem(42,-.1,1))
            fdom2 = DispersiveDomain(Interval(.1,.5),TwoLevelSystem(45,.3,2))
            lat = Cartesian(.001)
            latdom = LatticeDomain(bnd,lat)
            sim = Simulation(40,dom,latdom,fdom1,fdom2)
            nep = @test_nowarn HelmholtzNEP(sim)
            @test_nowarn nep(40)
            ω, ψ = @test_nowarn helmholtzeigen(nep,43;neigs=4,maxit=50)
        end
    end
end

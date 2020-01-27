using Iris
using NLsolve
using Test

@testset "Linear Scattering" begin
    @testset "1D" begin
        @testset "PML 2-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),PML,DirichletBC)
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            sim = Simulation(ω0, ndcavity, cell)
            ls = @test_nowarn HelmholtzLS(sim, ω0, 1)
            @test_nowarn scattering!(ls)
            @test ls.converged[]
            @test isapprox(ls.solution.scattered[288],-0.03032439146821593 + 0.09572749034950641im; atol=1e-3)
        end
        @testset "Matching 2-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),MatchedBC{1}(out=1),MatchedBC{2}(out=1))
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            sim = Simulation(ω0, ndcavity, cell)
            ls = @test_nowarn HelmholtzLS(sim, ω0, 1)
            @test_nowarn scattering!(ls)
            @test ls.converged[]
            @test isapprox(ls.solution.scattered[288],-0.03032439146821593 + 0.09572749034950641im; atol=1e-3)
        end
        @testset "PML 1-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),PML{1}(),DirichletBC)
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            sim = Simulation(ω0, ndcavity, cell)
            ls = @test_nowarn HelmholtzLS(sim, ω0, 1)
            @test_nowarn scattering!(ls)
            @test ls.converged[]
            @test isapprox(ls.solution.scattered[288],-0.1328681410970319 + 0.2692452251208401im; atol=1e-3)
        end
        @testset "Matched 1-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),PML{1}(),MatchedBC{1}(out=1),DirichletBC{2}())
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            sim = Simulation(ω0, ndcavity, cell)
            ls = @test_nowarn HelmholtzLS(sim, ω0, 1)
            @test_nowarn scattering!(ls)
            @test ls.converged[]
            @test isapprox(ls.solution.scattered[288],-0.1328681410970319 + 0.2692452251208401im; atol=1e-3)
        end
    end
end

@testset "Nonlinear Scattering" begin
    @testset "1D" begin
        @testset "PML 2-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),PML,DirichletBC)
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            dcavity = DispersiveDomain(Interval(-.3,.3),TwoLevelSystem(ω0,.1,2))
            sim = Simulation(ω0, ndcavity, cell, dcavity)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn fixedpoint(nls)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn nlsolve(nls)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test isapprox(nls.solution.scattered[288],-0.02742557704585827 + 0.1671878613596983im; atol=1e-3)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 10)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test isapprox(nls.solution.scattered[288],-0.24550578104277795 + 1.183068118874636im; atol=1e-3)
        end
        @testset "Matching 2-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),MatchedBC{1}(out=1),MatchedBC{2}(out=1))
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            dcavity = DispersiveDomain(Interval(-.3,.3),TwoLevelSystem(ω0,.1,2))
            sim = Simulation(ω0, ndcavity, cell, dcavity)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn fixedpoint(nls)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn nlsolve(nls)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test isapprox(nls.solution.scattered[288],-0.02742557704585827 + 0.1671878613596983im; atol=1e-3)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 10)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test isapprox(nls.solution.scattered[288],-0.24550578104277795 + 1.183068118874636im; atol=1e-3)
        end
        @testset "PML 1-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),PML{1}(),DirichletBC)
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            dcavity = DispersiveDomain(Interval(-.3,.3),TwoLevelSystem(ω0,.1,2))
            sim = Simulation(ω0, ndcavity, cell, dcavity)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn fixedpoint(nls)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn nlsolve(nls)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test isapprox(nls.solution.scattered[288],-0.22406937214601338 + 0.3040022689625704im; atol=1e-3)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 10)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test isapprox(nls.solution.scattered[288],-1.5947220620755382 + 2.8945059200579073im; atol=1e-3)
        end
        @testset "Matched 1-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),PML{1}(),MatchedBC{1}(out=1),DirichletBC{2}())
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            dcavity = DispersiveDomain(Interval(-.3,.3),TwoLevelSystem(ω0,.1,2))
            sim = Simulation(ω0, ndcavity, cell, dcavity)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn fixedpoint(nls)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn nlsolve(nls)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test isapprox(nls.solution.scattered[288],-0.22406937214601338 + 0.3040022689625704im; atol=1e-3)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 10)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test isapprox(nls.solution.scattered[288],-1.5947220620755382 + 2.8945059200579073im; atol=1e-3)
        end
    end
end

using Iris
using NLsolve
using Test
using LinearAlgebra

@testset "Linear Scattering" begin
    @testset "1D" begin
        @testset "PML 2-sided" begin
            ω=40
            sim = Iris.Examples.slab_pml(ω=ω, dx=.00035)
            ls = @test_nowarn HelmholtzLS(sim, ω, 1)
            @test_nowarn @inferred HelmholtzLS(sim, ω, 1)
            f(x) = x.operator
            @test_nowarn @inferred f(ls)
            f(x) = x.equivalent_source
            @test_nowarn @inferred f(ls)
            f(x) = x.j
            @test_nowarn @inferred f(ls)
            f(x) = x.solved
            @test_nowarn @inferred f(ls)
            f(x) = x.converged
            @test_nowarn @inferred f(ls)
            f(x) = x.solution
            @test_nowarn @inferred f(ls)
            f(x) = x.sim
            @test_nowarn @inferred f(ls)
            f(x) = x.ω
            @test_nowarn @inferred f(ls)
            f(x) = x.a
            @test_nowarn @inferred f(ls)
            f(x) = x.H
            @test_nowarn @inferred f(ls)
            f(x) = x.src
            @test_nowarn @inferred f(ls)
            f(x) = x.sol
            @test_nowarn @inferred f(ls)
            es = ls.equivalent_source
            f(es) = es.incoming_mask
            @test_nowarn @inferred f(es)
            f(es) = es.outgoing_mask
            @test_nowarn @inferred f(es)
            f(es) = es.field
            @test_nowarn @inferred f(es)
            f(es) = es.simulation
            @test_nowarn @inferred f(es)
            f(es) = es.ω
            @test_nowarn @inferred f(es)
            f(es) = es.a
            @test_nowarn @inferred f(es)
            f(es) = es.channelflux
            @test_nowarn @inferred f(es)
            sol = ls.solution
            f(sol) = sol.total
            @test_nowarn @inferred f(sol)
            f(sol) = sol.incident
            @test_nowarn @inferred f(sol)
            f(sol) = sol.scattered
            @test_nowarn @inferred f(sol)
            f(sol) = sol.ω
            @test_nowarn @inferred f(sol)
            f(sol) = sol.a
            @test_nowarn @inferred f(sol)
            @test_nowarn @inferred HelmholtzLS(sim, ω, 1)
            @test_nowarn scattering!(ls)
            @test ls.converged[]
            @test ls.helmholtz.A*ls.solution.total≈ls.j
            @test isapprox(ls.solution.scattered[288],0.07314178272991415 + 0.06502554346709624im; atol=1e-3)
        end
        @testset "Matching 2-sided" begin
            ω=40
            sim = Iris.Examples.slab_out(ω=ω, dx=.00035)
            ls = @test_nowarn HelmholtzLS(sim, ω, 1)
            @test_nowarn scattering!(ls)
            @test ls.converged[]
            @test ls.helmholtz.A*ls.solution.total≈ls.j
            @test isapprox(ls.solution.scattered[288],0.07314178272991415 + 0.06502554346709624im; atol=1e-3)
        end
        @testset "PML 1-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),PML{1}(),DirichletBC)
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            sim = Simulation(ω0, ndcavity, cell)
            smooth!(sim)
            ls = @test_nowarn HelmholtzLS(sim, ω0, 1)
            @test_nowarn scattering!(ls)
            @test ls.converged[]
            @test ls.helmholtz.A*ls.solution.total≈ls.j
            @test isapprox(ls.solution.scattered[288],-0.1290389888137668 + 0.27289178643785555im; atol=1e-3)
        end
        @testset "Matched 1-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),PML{1}(),MatchedBC{1}(out=1),DirichletBC{2}())
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            sim = Simulation(ω0, ndcavity, cell)
            smooth!(sim)
            ls = @test_nowarn HelmholtzLS(sim, ω0, 1)
            @test_nowarn scattering!(ls)
            @test ls.converged[]
            @test ls.helmholtz.A*ls.solution.total≈ls.j
            @test isapprox(ls.solution.scattered[288],-0.1290389888137668 + 0.27289178643785555im; atol=1e-3)
        end
    end
end

@testset "Nonlinear Scattering" begin
    @testset "1D" begin
        @testset "PML 2-sided" begin
            ω=40
            sim = Iris.Examples.slab_pml(ω=ω, dx=.00035, D=.1, γ=2)
            nls = @test_nowarn HelmholtzNLS(sim, ω, 1)
            @test_nowarn @inferred HelmholtzNLS(sim,ω,1)
            f(x) = x.linearscatter
            @test_nowarn @inferred f(nls)
            f(x) = x.ψ
            @test_nowarn @inferred f(nls)
            f(x) = x.residual
            @test_nowarn @inferred f(nls)
            f(x) = x.x
            @test_nowarn @inferred f(nls)
            f(x) = x.fixedpoint
            @test_nowarn @inferred f(nls)
            f(x) = x.lupack
            @test_nowarn @inferred f(nls)
            f(x) = x.a
            @test_nowarn @inferred f(nls)
            f(x) = x.equivalent_source
            @test_nowarn @inferred f(nls)
            f(x) = x.helmholtz
            @test_nowarn @inferred f(nls)
            f(x) = x.operator
            @test_nowarn @inferred f(nls)
            f(x) = x.j
            @test_nowarn @inferred f(nls)
            f(x) = x.simulation
            @test_nowarn @inferred f(nls)
            f(x) = x.solution
            @test_nowarn @inferred f(nls)
            f(x) = x.ω
            @test_nowarn @inferred f(nls)
            f(x) = x.results
            @test_nowarn scattering!(nls)
            @test_nowarn @inferred f(nls)
            nls = @test_nowarn HelmholtzNLS(sim, ω, 1)
            @test_nowarn fixedpoint(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            nls = @test_nowarn HelmholtzNLS(sim, ω, 1)
            @test_nowarn nlsolve(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            nls = @test_nowarn HelmholtzNLS(sim, ω, 1)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            @test isapprox(nls.solution.scattered[288],0.1881169235362705 + 0.15566101877014im; atol=1e-3)
            nls = @test_nowarn HelmholtzNLS(sim, ω, 10)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            @test isapprox(nls.solution.scattered[288], 1.1192308728340543 + 0.7280855538638448im; atol=1e-3)
        end
        @testset "Matching 2-sided" begin
            ω0=40
            sim = Iris.Examples.slab_out(ω=ω0, dx=.00035, D=.1, γ=2)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn fixedpoint(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn nlsolve(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            @test isapprox(nls.solution.scattered[288],0.18811652728844638 + 0.15566252089682026im; atol=1e-3)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 10)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            @test isapprox(nls.solution.scattered[288],1.1192537618187042 + 0.7280915535384698im; atol=1e-3)
        end
        @testset "PML 1-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-1,1),PML{1}(),MatchedBC{1}(out=1),DirichletBC{2}())
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            dcavity = DispersiveDomain(Interval(-.3,.3),TwoLevelSystem(ω0,.1,2))
            sim = Simulation(ω0, ndcavity, cell, dcavity)
            smooth!(sim)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn fixedpoint(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn nlsolve(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            @test isapprox(nls.solution.scattered[288],-0.07570911698903976 + 0.3176555836644436im; atol=1e-3)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 10)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            @test isapprox(nls.solution.scattered[288],0.5081023860501906 + 2.6581334068616136im; atol=1e-3)
        end
        @testset "Matched 1-sided" begin
            ω0 = 40
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-1,1),PML{1}(),MatchedBC{1}(out=1),DirichletBC{2}())
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            dcavity = DispersiveDomain(Interval(-.3,.3),TwoLevelSystem(ω0,.1,2))
            sim = Simulation(ω0, ndcavity, cell, dcavity)
            smooth!(sim)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn fixedpoint(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn nlsolve(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 1)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            @test isapprox(nls.solution.scattered[288],-0.07570911698903976 + 0.3176555836644436im; atol=1e-3)
            nls = @test_nowarn HelmholtzNLS(sim, ω0, 10)
            @test_nowarn scattering!(nls)
            @test nls.converged[]
            @test nls.helmholtz.A*nls.solution.total≈nls.j
            @test isapprox(nls.solution.scattered[288],0.5081023860501906 + 2.6581334068616136im; atol=1e-3)
        end
    end

    @testset "SPA" begin
        ω=100.4
        sim = Iris.Examples.slab_pml(ω=ω, dx=.00035, D=-1, γ=2)
        nls = @test_nowarn HelmholtzNLS(sim, ω, 50, 10)
        scattering!(nls)
        @test_nowarn @inferred spa(nls; nev=2)
    end
end

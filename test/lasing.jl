using Iris
using NLsolve
using Test

@testset "SALT" begin
    @testset "1D PML" begin
        ω0 = 150.5
        nD = 153
        D0s = LinRange(0,.9,nD)
        global flag = false
        global Ω = ω0
        for D ∈ D0s
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),PML,DirichletBC)
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            tls = TwoLevelSystem(ωa=ω0,D0=D,γp=1)
            dcavity = DispersiveDomain(Interval(-.5,.5),tls)
            sim = Simulation(ω0, ndcavity, dcavity, cell)
            if !flag
                nep = HelmholtzNEP(sim)
                ωs, ψs = helmholtzeigen(nep, Ω; neigs=1, maxit=40)
                global Ω = ωs[1]
                if imag(ωs[1]) > 0
                    global flag = true
                end
            end
            if flag
                salt = @test_nowarn HelmholtzSALT(sim, 1)
                @test_nowarn nlsolve(salt, real(Ω), ψs; show_trace=false, ftol=1e-6, iterations=250)
                global ψs = salt.ψs
                @test salt.converged[]
            end
        end
    end
    @testset "1D Matching" begin
        ω0 = 150.5
        nD = 153
        D0s = LinRange(0,.9,nD)
        global flag = false
        global Ω = ω0
        for D ∈ D0s
            lat = Lattices.Cartesian(.00035)
            bnd = Boundary(Interval(-.7,.7),MatchedBC{1}(out=1),MatchedBC{2}(out=1))
            cell = LatticeDomain(bnd,lat)
            ndcavity = NondispersiveDomain(Interval(-.5,.5),3)
            tls = TwoLevelSystem(ωa=ω0,D0=D,γp=1)
            dcavity = DispersiveDomain(Interval(-.5,.5),tls)
            sim = Simulation(ω0, ndcavity, dcavity, cell)
            if !flag
                nep = HelmholtzNEP(sim)
                ωs, ψs = helmholtzeigen(nep, Ω; neigs=1, maxit=100)
                global Ω = ωs[1]
                if imag(ωs[1]) > 0
                    global flag = true
                end
            end
            if flag
                salt = @test_nowarn HelmholtzSALT(sim, 1)
                @test_nowarn nlsolve(salt, real(Ω), ψs; show_trace=false, ftol=1e-6, iterations=250)
                global ψs = salt.ψs
                @test salt.converged[]
            end
        end
    end
end

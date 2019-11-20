## LEP
using Iris
using Plots
lat = Lattice(2π/1564)
    bnd1 = Boundary(Interval(0,3),noBL,DirichletBC;depth=.6)
    bnd2 = Boundary(Interval(1,2))
    dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0);name=:air)
    dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1),TwoLevelSystem(30,-.1,1);name=:cavity,type=:dielectric)
    sim = Simulation(30,dom2,dom1)
lep = MaxwellLEP(sim)
ωl, ψl = maxwelleigen(lep,30;nev=3,lupack=PSolver())
display(ωl); plot(ψl)

mb = Iris.TimeDomain.maxwell_bloch(sim,30)
    # mb.fields.E.y[12:15] .= 1
    # mb.fields.B.z[12:15] .= 1
    mb.fields.E.y[50 .+ (-40:40)] = exp.(-((-40:40)./10).^2)
    mb.fields.B.z[50 .+ (-40:40)] = exp.(-((-40:40)./10).^2)
    plot(mb.fields.E)
    plot!(mb.fields.B)
mb(210000)
    plot(mb.fields.E,ylims=(0,1))
plot(mb.fields.B,ylims=(0,1))

## CF
using Iris
using Plots
lat = Lattice(π/30/15/3)
bnd1 = Boundary(Interval(0,3),PML,DirichletBC;depth=.5)
bnd2 = Boundary(Interval(1,2))
dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1),TwoLevelSystem(30,-.1,1))
sim = Simulation(30,dom2,dom1)
cf = MaxwellCF(sim)
η, u = maxwelleigen(cf,30;nev=2,verbose=true)
display(η); plot(u)

## NEP
using Iris
using NonlinearEigenproblems
using Plots
lat = Lattice(.005)
bnd1 = Boundary(Interval(0.5,2.5),MatchedBC{1}(in=[1]),MatchedBC{2}(in=[1]))
# bnd1 = Boundary(Interval(0,3),PML,DirichletBC;depth=.5)
bnd2 = Boundary(Interval(1,2))
dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1),TwoLevelSystem(30,-.15,1))
sim = Simulation(30,dom2,dom1)
nep = maxwell_nep(sim)
ωnl, ψnl = maxwell_eigen(nep,30;neigs=4,maxit=51)
display(ωnl); plot(ψnl)

## SALT
using Iris
using Plots
using NLsolve
# lat = Lattice(2π/30/15/30)
# bnd1 = Boundary(Interval(0,3),MatchedBC{2}(out=[1]),MatchedBC{1}(out=[1]))
bnd1 = Boundary(Interval(0,3),PML,DirichletBC;depth=.8)
bnd2 = Boundary(Interval(1,2))
dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1),TwoLevelSystem(15,.15,1))
sim = Simulation(15,dom2,dom1)
salt = maxwell_salt(sim,1)
nlsolve(salt,real(ωnl[2]),ψnl(2);show_trace=true)
plot(salt)

## Linear Scattering
using Iris
using Plots
lat = Lattice(2π/30/15/4)
bnd1 = Boundary(Interval(0,3),MatchedBC{1}(out=[1]),MatchedBC{2}(out=[1]))
bnd2 = Boundary(Interval(1,2))
dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
dom2 = Domain(bnd2,lat,DielectricFunction(4),PumpFunction(1),TwoLevelSystem(30.17,-.35,.12))
sim = Simulation(30,dom2,dom1)
@time sct = scattering(sim,29.6,0,1)
plot(sct)

S1 = smatrix(sim,30.2.+collect(LinRange(-2,2,51)))
    plot(S1)
S2 = smatrix(sim,30.2.+LinRange(-2,2,51))
    plot(S2)

## Nonlinear Scattering
using Iris
using Plots
using NLsolve
lat = Lattice(2π/30/15/3)
bnd1 = Boundary(Interval(.5,2.5),MatchedBC{2}(out=[1]),MatchedBC{1}(out=[1]))
bnd2 = Boundary(Interval(1,2))
dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1),TwoLevelSystem(30,.14,1))
sim = Simulation(30,dom2,dom1)
a=2.243; θ=π
    nls = maxwell_nls(sim,scpa.ωs[1],a,cis(θ)*a)
    nlsolve(nls;show_trace=true)
    plot(nls.sol.total,ylims=(0,.19))
plot!(scpa)


## CPA
using Iris
using Plots
using NLsolve
lat = Lattice(2π/10/2)
bnd1 = Boundary(Interval(.5,2.5),MatchedBC{2}(in=[1]),MatchedBC{1}(in=[1]))
bnd2 = Boundary(Interval(1,2))
dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1),TwoLevelSystem(30,-.35,1))
sim = Simulation(30,dom2,dom1)
@time scpa = maxwell_scpa(sim,1)
nlsolve(scpa,real(ωnl[1]),ψnl(1);show_trace=true)
plot(scpa)

## NLSCATTER
lat = Lattice(2π/30/15/3)
bnd1 = Boundary(Interval(0,3),MatchedBC{2}(out=[1]),MatchedBC{1}(out=[1]))
bnd2 = Boundary(Interval(1,2))
dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1),TwoLevelSystem(30,-.15,1))
sim = Simulation(scpa.ωs[1],dom2,dom1)

N = 20
a = LinRange(.01,2,N)
sct = []
for i ∈ eachindex(a)
    if i==1
        push!(sct,scattering_nl(sim,scpa.ωs[1],a[1]*scpa.ψs.left[2],a[1]*scpa.ψs.right[2]))
    else
        push!(sct,scattering_nl(sim,scpa.ωs[1],a[i]*scpa.ψs.left[2],a[1]*scpa.ψs.right[2];ψ_init=sct[i-1].tot))
    end
end
_,indL = findmin(map(s->s.x,sim.x))
_,indR = findmax(map(s->s.x,sim.x))
RL = Vector{ComplexF64}(undef,length(a))
RR = similar(RL)
# T = similar(RL)
for i ∈ eachindex(RL)
    RL[i] = sct[i].sct.y[indL]/sqrt(scpa.ωs[1])
    RR[i] = sct[i].sct.y[indR]/sqrt(scpa.ωs[1])
    # T[i] = sct[i].tot.y[indR]/sqrt(scpa.ωs[1])
end
plot(2abs2.(a),abs2.(RL)+abs2.(RR))

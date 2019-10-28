## test Point
using Iris
p = Point(1)
p = Point(1,2)
p = Point((1,2,3))

p.x
p.y
p.ϕ
p.r

@code_warntype Point(1)
@code_warntype Point(1,2)
@code_warntype Point((1,2,3))
f(p) = p.x
@code_warntype f(p)


## test BoundaryLayers
using Iris
PML{1}(.3)
cPML{2}(.3)
noBL{1}(.4)
conj(PML{3}(.3))

pml = PML{1}(.3,0,.3)
@code_warntype pml(.1)

nbl = noBL{1}(.3,0,.3)
@code_warntype nbl(.1)


## test Shapes1D
using Iris
interval = Interval(0,.3)
@code_warntype interval(.2)


## test DielectricFunction
using Iris
DielectricFunction()
de = DielectricFunction(1,.1)
desin(p::Point,params::Dict) = complex(params[:n1],params[:n2]*sin(params[:period]*p.x))
DielectricFunction(desin,2,.3)
de2 = DielectricFunction(desin,Dict(:period=>π,:n1=>2,:n2=>.1))
de2(Point(2.1))
F = PumpFunction()
F(Point(1.4))
fsin(p::Point,params::Dict) = params[:f0] + params[:f1]*sin(params[:period]*p.x)
F2 = PumpFunction(fsin,Dict(:period=>π,:f0=>1,:f1=>.1))
F2(Point(1.1))

@code_warntype DielectricFunction()
@code_warntype DielectricFunction(1,.1)
@code_warntype DielectricFunction(desin,2,.3)
@code_warntype DielectricFunction(desin,Dict{Symbol,Float64}(:period=>π,:n1=>2,:n2=>.1))
@code_warntype de2(Point(2.1))
@code_warntype PumpFunction()
@code_warntype F(Point(1.4))
@code_warntype PumpFunction(fsin,Dict(:period=>π,:f0=>1,:f1=>.1))
@code_warntype F2(1.1)


## lattices
using Iris
lat = Lattice(.01)
Lattice(1.0,1.0;type=:Polar)
lat = Lattice(1.0,1.0,1.0;angles=(.1,.1,0.))

@code_warntype Lattice(1.0)
@code_warntype Lattice(1.0,1.0)
@code_warntype Lattice(1.0,1.0,1.0;angles=(.1,.1,0.))
f(lat) = lat.r0
@code_warntype f(lat)


## Boundary
using Iris
interval = Interval(0,1)

Boundary(interval,NeumannBC,PML)
@code_warntype Boundary(interval,NeumannBC,PML)
Boundary(interval,NeumannBC)
@code_warntype Boundary(interval,NeumannBC)
Boundary(interval;depth=.2)
@code_warntype Boundary(interval;depth=.2)

Boundary(NeumannBC,interval,PML)
@code_warntype Boundary(NeumannBC,interval,PML)
Boundary(NeumannBC,interval)
@code_warntype Boundary(NeumannBC,interval)
Boundary(PML,interval)
@code_warntype Boundary(PML,interval)

Boundary(NeumannBC,PML,interval)
@code_warntype Boundary(NeumannBC,PML,interval)
Boundary(PML,NeumannBC,interval)
@code_warntype Boundary(PML,NeumannBC,interval)

bnd = Boundary(interval,DirichletBC{2}())
@code_warntype Boundary(interval,DirichletBC{2}())

bcs = (DirichletBC{1}(),DirichletBC{2}())
bls = (PML{1}(.3),PML{2}(.3))
bnd = Boundary(interval,bcs,bls)

@code_warntype Boundary(interval,bcs,bls)


## Domain
using Iris
lat = Lattice(.013)
interval = Interval(0,1)
bcs = (DirichletBC{1}(),DirichletBC{2}())
bls = (PML{1}(.3),PML{2}(.3))
bnd = Boundary(interval,bcs,bls)

dom = Domain(bnd,lat,DielectricFunction(1),PumpFunction(0))
@code_warntype Domain(bnd,lat,DielectricFunction(1),PumpFunction(0))


## Simulation
using Iris
lat = Lattice(.013)
interval1 = Interval(0,2)
interval2 = Interval(.5,1)
bcs = (DirichletBC{1}(),DirichletBC{2}())
bls = (PML{1}(.2),PML{2}(.2))
bnd1 = Boundary(interval1,bcs,bls)
bnd2 = Boundary(interval2,bcs,bls)

dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
dom2 = Domain(bnd2,lat,DielectricFunction(2),PumpFunction(0))
sim = Simulation(10,dom1,dom2)


##
using NLsolve
using Iris
    using Plots
    lat = Lattice(2π/30/15/3)
    interval1 = Interval(0,3)
    interval2 = Interval(1,2)
    bcs1 = (MatchedBC{1}(out=[1]),MatchedBC{2}(out=[1]))
    bls1 = (noBL{1}(4*2π/30),noBL{2}(4*2π/30))
    bcs2 = (noBC{1}(),noBC{2}())
    bls2 = (noBL{1}(),noBL{2}())
    bnd1 = Boundary(interval1,bcs1,bls1)
    bnd2 = Boundary(interval2,bcs2,bls2)

    dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
    dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1),TwoLevelSystem(30,-.15,2))
    sim = Simulation(10,dom2,dom1)

    m = maxwell(sim)

    @time ωl,ψl = maxwell_eigs_l(sim,30;nev=3)
    sct = scattering_nl(sim,30;show_trace=true)
plot(sctl)
plot(sct)

@time ms = maxwell_salt(sim,1)
@time mc = maxwell_scpa(sim,1)

@time nlsolve(ms,real(ωl),ψl;show_trace=true);
@time nlsolve(mc,real(ωl),.1ψl;show_trace=true);
@time SALT(sim,real(ωl),ψl;show_trace=true)

@time sct = scattering(sim,10;start=.5,stop=2.5)

@btime eig_kl(sim,30);
@code_warntype eig_kl(sim,30)
using NonlinearEigenproblems
nep = maxwell_nep(sim)
ω,ψ = contour_beyn(nep;σ=8,radius=.3,N=75,neigs=2,k=3); println(ω)
@btime eig_knl(sim,30; maxit=40);
@code_warntype eig_knl(sim,30; maxit=40)

η,u = eig_cf(sim,30)
@btime eig_cf(sim,30);
@code_warntype eig_cf(sim,30)

perm = sortperm(map(a->a.vec[1],sim.x))
plot(map(a->a.vec[1],sim.x)[perm],100abs2.(ψl[end÷3+1:2end÷3,1])[perm])
plot(map(a->a.vec[1],sim.x)[perm],100abs2.(ψnl[end÷3+1:2end÷3,1])[perm])
plot(map(a->a.vec[1],sim.x)[perm],100abs2.(u[end÷3+1:2end÷3,1])[perm])


## dispersions
using Iris
using NonlinearEigenproblems
lat = Lattice(2π/30/30/3)
interval1 = Interval(0,3)
interval2 = Interval(1,2)
bcs1 = (MatchedBC{1}(out=[1]),MatchedBC{2}(out=[1]))
bls1 = (noBL{1}(4*2π/30),noBL{2}(4*2π/30))
bcs2 = (noBC{1}(),noBC{2}())
bls2 = (noBL{1}(),noBL{2}())
bnd1 = Boundary(interval1,bcs1,bls1)
bnd2 = Boundary(interval2,bcs2,bls2)

dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1),TwoLevelSystem(31,.1,2))
sim = Simulation(30,dom2,dom1)

nep = maxwell_nep(sim)
ω,ψ = contour_beyn(nep;σ=32,radius=1,neigs=3,k=5,logger=2,N=75)
ω2,ψ2 = contour_beyn(nep2;σ=32,radius=1,neigs=3,k=5,logger=2,N=75)
# newton(nep1;λ=32,maxit=50,logger=2)
# newton(nep2;λ=32,maxit=50,logger=2)
# augnewton(nep1;λ=32,maxit=50,logger=2)
# augnewton(nep2;λ=32,maxit=50,logger=2)
resinv(nep.nep;λ=31.3,logger=2,maxit=200)
resinv(nep.nep;λ=32,logger=2)
quasinewton(nep.nep;λ=32,logger=2)
quasinewton(nep.nep;λ=32,logger=2,maxit=200)
mslp(nep.nep;λ=32,logger=2,linsolvercreator=IrisLinSolverCreator(USolver))
mslp(nep.nep;λ=32,logger=2)
# implicitdet(nep1;λ=32,maxit=150,logger=2)
# implicitdet(nep2;λ=32,maxit=150,logger=2)
# nlar(nep1;λ=32,logger=2,neigs=2)
# nlar(nep2;λ=32,logger=2,neigs=2)
# jd_betcke(nep1;λ=32,logger=2)
# jd_betcke(nep2;λ=32,logger=2)
# jd_effenberger(nep2;λ=32,logger=2)
# jd_effenberger(nep2;λ=32,logger=2)
# contour_block_SS(nep1;σ=32,neigs=1,radius=1,k=6,logger=2,N=150)
# contour_block_SS(nep2;σ=32,neigs=3,k=5,logger=2,N=75)
a,b = iar(nep.nep;σ=32,logger=1,neigs=2)
iar(nep2;σ=32,logger=1,neigs=2)
iar_chebyshev(nep.nep;σ=32,logger=1,neigs=2)
iar_chebyshev(nep2;σ=32,logger=1,neigs=2)
tiar(nep1;σ=32,logger=1,neigs=2)
tiar(nep2;σ=32,logger=1,neigs=2)
# ilan(nep1;σ=32,logger=1,neigs=1)
# ilan(nep2;σ=32,logger=1,neigs=1)


## test jacobian

using LinearAlgebra
using SparseArrays

tls = TwoLevelSystem(10,1,1)
ωs = [9. 11.]

dω = .00001
v = [1 5 .1]
Dω = dω*v

γ1 = Iris.Common.Dispersions.TwoLevelSystems.γ(tls,ωs,ψl)
γ2 = Iris.Common.Dispersions.TwoLevelSystems.γ(tls,ωs+Dω,ψl)
dγ = (γ2-γ1)/dω
∂γ = Iris.Common.Dispersions.TwoLevelSystems.∂γ∂ω(tls,ωs,ψl).*v
norm(dγ-∂γ)

Γ1 = Iris.Common.Dispersions.TwoLevelSystems.Γ(tls,ωs,ψl)
Γ1-abs2.(γ1)
Γ2 = Iris.Common.Dispersions.TwoLevelSystems.Γ(tls,ωs+Dω,ψl)
dΓ = (Γ2-Γ1)/dω
∂Γ = Iris.Common.Dispersions.TwoLevelSystems.∂Γ∂ω(tls,ωs,ψl).*v
norm(dΓ-∂Γ)

h1 = Iris.Common.Dispersions.TwoLevelSystems.H(tls,ωs,ψl)
h2 = Iris.Common.Dispersions.TwoLevelSystems.H(tls,ωs+Dω,ψl)
dh = (h2 - h1)/dω
∂h = Iris.Common.Dispersions.TwoLevelSystems.∂H∂ω(tls,ωs,ψl)
norm(∂h[:,1]*v[1]+∂h[:,2]*v[2]+∂h[:,3]*v[3] - dh)

χ1 = Iris.Common.Dispersions.TwoLevelSystems.χ(tls,ωs,ψl)
χ2 = Iris.Common.Dispersions.TwoLevelSystems.χ(tls,ωs+Dω,ψl)
dχψ = ((χ2-χ1)/dω)*ψl[:]
∂χψ = Iris.Common.Dispersions.TwoLevelSystems.∂χ∂ωψ(tls,ωs,ψl)
norm(∂χψ[1][:]*v[1] + ∂χψ[2][:]*v[2] + ∂χψ[3][:]*v[3] - dχψ)


####

using Iris

lat = Lattice(2π/30/10/3)
    interval1 = Interval(0,3)
    interval2 = Interval(1,2)
    bcs1 = (MatchedBC{1}(out=[1]),MatchedBC{2}(out=[1]))
    bls1 = (noBL{1}(4*2π/30),noBL{2}(4*2π/30))
    bcs2 = (noBC{1}(),noBC{2}())
    bls2 = (noBL{1}(),noBL{2}())
    bnd1 = Boundary(interval1,bcs1,bls1)
    bnd2 = Boundary(interval2,bcs2,bls2)

    dom1 = Domain(bnd1,lat,DielectricFunction(1),PumpFunction(0))
    dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1))
    sim = Simulation(30,dom2,dom1)

using Random
ωs = [8. 9. 11.]
ψ = rand(MersenneTwister(1),ComplexF64,3length(sim.x),length(ωs))
d = .0001
dψ = d*rand(ComplexF64,3length(sim.x),length(ωs))
tls = TwoLevelSystem(10,.2,1)
χ1 = susceptability(tls,ωs,ψ)
χ2 = susceptability(tls,ωs,ψ+dψ)
dχ = (χ2-χ1)
J = jacobian_lasing(tls,ωs,ψ)
J*vcat(real(dψ[:]),imag(dψ[:]))

using Iris


## test Point

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

PML{1}(.3)
cPML{2}(.3)
noBL{1}(.4)
conj(PML{3}(.3))

pml = PML{1}(.3,0,.3)
@code_warntype pml(.1)

nbl = noBL{1}(.3,0,.3)
@code_warntype nbl(.1)

## test Shapes1D

interval = Interval(0,.3)
@code_warntype interval(.2)

## test DielectricFunction

DielectricFunction()
de = DielectricFunction(1,.1)
desin(x::Number,params::Dict) = complex(params[:n1],params[:n2]*sin(params[:period]*x))
DielectricFunction(desin,2,.3)
de2 = DielectricFunction(desin,Dict(:period=>π,:n1=>2,:n2=>.1))
de2(2.1)
F = PumpFunction()
F(1.4)
fsin(x::Number,params::Dict) = params[:f0] + params[:f1]*sin(params[:period]*x)
F2 = PumpFunction(fsin,Dict(:period=>π,:f0=>1,:f1=>.1))
F2(1.1)

@code_warntype DielectricFunction()
@code_warntype DielectricFunction(1,.1)
@code_warntype DielectricFunction(desin,2,.3)
@code_warntype DielectricFunction(desin,Dict{Symbol,Float64}(:period=>π,:n1=>2,:n2=>.1))
@code_warntype de2(2.1)
@code_warntype PumpFunction()
@code_warntype F(1.4)
@code_warntype PumpFunction(fsin,Dict(:period=>π,:f0=>1,:f1=>.1))
@code_warntype F2(1.1)

## lattices

lat = Lattice(.01)
# Lattice(1.0,1.0;type=:Polar)
# lat = Lattice(1.0,1.0,1.0;angles=(.1,.1,0.))

@code_warntype Lattice(1.0)
@code_warntype Lattice(1.0,1.0)
@code_warntype Lattice(1.0,1.0,1.0;angles=(.1,.1,0.))
f(lat) = lat.r0
@code_warntype f(lat)

## Boundary

interval = Interval(0,1)
bcs = (DirichletBC{1}(),DirichletBC{2}())
bls = (PML{1}(.3),PML{2}(.3))
bnd = Boundary(interval,bcs,bls)

@code_warntype Boundary(interval,bcs,bls)

##

lat = Lattice(.013)
interval = Interval(0,1)
bcs = (DirichletBC{1}(),DirichletBC{2}())
bls = (PML{1}(.3),PML{2}(.3))
bnd = Boundary(interval,bcs,bls)

dom = Domain(bnd,lat,DielectricFunction(1),PumpFunction(0))
@code_warntype Domain(bnd,lat,DielectricFunction(1),PumpFunction(0))

## Simulation

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
using Iris
using NonlinearEigenproblems
using Plots
# using Pardiso

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
dom2 = Domain(bnd2,lat,DielectricFunction(3),PumpFunction(1))
sim = Simulation(30,dom2,dom1)

ωl,ψl = eig_kl(sim,30)
@btime eig_kl(sim,30);
@code_warntype eig_kl(sim,30)

ωnl,ψnl = eig_knl(sim,30; maxit=40)
@btime eig_knl(sim,30; maxit=40);
@code_warntype eig_knl(sim,30; maxit=40)

η,u = eig_cf(sim,30)
@btime eig_cf(sim,30);
@code_warntype eig_cf(sim,30)

perm = sortperm(map(a->a.vec[1],sim.x))
plot(map(a->a.vec[1],sim.x)[perm],100abs2.(ψl[end÷3+1:2end÷3,1])[perm])
plot(map(a->a.vec[1],sim.x)[perm],100abs2.(ψnl[end÷3+1:2end÷3,1])[perm])
plot(map(a->a.vec[1],sim.x)[perm],100abs2.(u[end÷3+1:2end÷3,1])[perm])

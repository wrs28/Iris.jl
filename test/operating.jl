## Simulation
using Iris
using Plots
lat = Lattice(.005)
bnd = Boundary(Interval(-1,1.1),PML)
dom = LatticeDomain(bnd,lat)
dom1 = NondispersiveDomain(Interval(-.5,.5),3)
tls = TwoLevelSystem(D0=.1,omega=30,gamma=1)
dom2 = DispersiveDomain(Interval(-.5,.5),tls)
sim = Simulation(30,dom,dom1,dom2)
smooth!(sim)
plot(sim)

## LEP
lep = MaxwellLEP(sim)
ωl, ψl = maxwelleigen(lep,30;nev=2)
display(ωl); plot(ψl)

## CF
cf = MaxwellCF(sim)
η, u = maxwelleigen(cf,30)
display(η); plot(u)

## Linear Scattering
ls1 = MaxwellLS(sim,30,1,0)
ls2 = MaxwellLS(sim,30,0,1)
scattering!(ls1)
scattering!(ls2)
plot(ls1)
plot(ls2)

## SMatrix
S = smatrix(sim,LinRange(29,31,201))
plot(S)

##
using NLsolve
nls1 = MaxwellNLS(sim,30,1,0)
fixedpoint(nls1; show_trace=true)

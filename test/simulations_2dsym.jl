##

using Iris
using Plots

shape = Square(2,0,0)
bnd = Boundary(shape,DirichletBC,PML)
lat = Lattices.Cartesian(.1,.125)
cell = LatticeDomain(bnd,lat)

cavity = NondispersiveDomain(Circle(.5,0,0),3)

@code_warntype Simulation(30,cell,cavity)
@btime Simulation(30,cell,cavity)
sim = Simulation(30,cell,cavity)

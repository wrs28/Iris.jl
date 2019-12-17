##

using Iris
using Plots

shape = Square(1,1.1,0)
plot(shape)

bnd = Boundary(shape,DirichletBC,PML)
plot(bnd)

lat = Lattices.Cartesian(.1,.125; angles=.1)
plot(lat)

@code_warntype LatticeDomain(bnd,lat)
@btime dom = LatticeDomain(bnd,lat)
dom = LatticeDomain(bnd,lat)
plot(LatticeDomain(bnd,lat))

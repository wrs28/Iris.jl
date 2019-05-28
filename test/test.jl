using Revise, Plots
using Iros

c = Circle(1,2,3)
plot(c)
e = Ellipse(1,2,3,4,.1)
plot(e)
d = DeformedDisk{2}(1,2,3,[4,5],[.06,.07],[8,9],.1)
plot(d)

a = Annulus(1,2,3,4,5)
plot(a)

s = Square(1,2,3,.1)
plot(s)

r = Rectangle(1,2,3,4,.1)
plot(r)

lat = Lattices.Cartesian(1/11,1/11)

DirichletBC()

bnd = Boundary(DeformedDisk{2}(1,1,2,[2,3],[.5,.1],[0,0]),DirichletBC())

dom = Cavity(bnd,lat)

sim = Simulation(dom)

# bnd = Boundary(Rectangle(1,.1,0,0,0.1),DirichletBC())
lat = Lattices.Cartesian(1/101,1/101,.1)
dom = Cavity(bnd,lat)
sim = Simulation(dom)
plot(sim)

k,ψ = eig_kl(sim,5,nev=25)


n=0
n+=1;print(n,": ");println(k[n]);plot(sim,ψ,n,by=abs2)
    # plot!(bnd)

plot(xlims=(0.4,0.5),ylims=(2.7,2.8))
    # plot!(sim,ψ,3,by=abs2,ms=5)
    scatter!(sim.x[sim.surface],sim.y[sim.surface],label="Surface")
    # scatter!(sim.x[sim.exterior],sim.y[sim.exterior],label="Exterior")
    indp = ind+1
    inds = findall(sim.link_y_bool.&(sim.link_start.==indp))
    scatter!(sim.x[[indp]],sim.y[[indp]],label="point",shape=:star,ms=10)
    scatter!(sim.x[sim.link_stop[inds]],sim.y[sim.link_stop[inds]],marker_z=sim.link_weight[inds]*lat.dx*lat.dy,color=:diverging)
    plot!(clims=(-1,1))
    plot!(bnd,label=:bnd)
    plot!(leg=true)



plot!(sim,ψ,6,by=abs2,ms=5)
    # scatter(sim.x[sim.interior],sim.y[sim.interior])
scatter!(sim.x[sim.surface],sim.y[sim.surface])
scatter!(sim.x[sim.exterior],sim.y[sim.exterior])
    # scatter!(sim.x[sim.bulk],sim.y[sim.bulk])
    scatter!(sim.x[sim.corner],sim.y[sim.corner])
    # plot!(sim.tessellation)
    # plot!(bnd.shape)

    scatter!(sim.x[sim.lattice_wall],sim.y[sim.lattice_wall])

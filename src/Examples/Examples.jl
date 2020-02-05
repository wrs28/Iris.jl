module Examples

using ..Common

function slab_pml(; D::Real=0, ω::Real = 40, γ::Real=2, n::Real=3, dx::Real=2π/ω/n/20)
    cavity = NondispersiveDomain(Interval(-.5,.5),n)
    pump = DispersiveDomain(Interval(-.5,.5),TwoLevelSystem(ω,D,γ))
    bnd = Boundary(Interval(-1,1),PML,DirichletBC)
    lat = Cartesian(dx)
    latdom = LatticeDomain(bnd,lat)
    return smooth!(Simulation(ω,cavity,pump,latdom))
end

function slab_cpml(; D::Real=0, ω::Real = 40, γ::Real=2, n::Real=3, dx::Real=2π/ω/n/20)
    cavity = NondispersiveDomain(Interval(-.5,.5),n)
    pump = DispersiveDomain(Interval(-.5,.5),TwoLevelSystem(ω,D,γ))
    bnd = Boundary(Interval(-1,1),cPML,DirichletBC)
    lat = Cartesian(dx)
    latdom = LatticeDomain(bnd,lat)
    return smooth!(Simulation(ω,cavity,pump,latdom))
end

function slab_out(; D::Real=0, ω::Real = 40, γ::Real=2, n::Real=3, dx::Real=2π/ω/n/20)
    cavity = NondispersiveDomain(Interval(-.5,.5),n)
    pump = DispersiveDomain(Interval(-.5,.5),TwoLevelSystem(ω,D,γ))
    bnd = Boundary(Interval(-1,1),MatchedBC{1}(out=[1]),MatchedBC{2}(out=[1]))
    lat = Cartesian(dx)
    latdom = LatticeDomain(bnd,lat)
    return smooth!(Simulation(ω,cavity,pump,latdom))
end

function slab_in(; D::Real=0, ω::Real = 40, γ::Real=2, n::Real=3, dx::Real=2π/ω/n/20)
    cavity = NondispersiveDomain(Interval(-.5,.5),n)
    pump = DispersiveDomain(Interval(-.5,.5),TwoLevelSystem(ω,D,γ))
    bnd = Boundary(Interval(-1,1),MatchedBC{1}(in=[1]),MatchedBC{2}(in=[1]))
    lat = Cartesian(dx)
    latdom = LatticeDomain(bnd,lat)
    return smooth!(Simulation(ω,cavity,pump,latdom))
end

end # module

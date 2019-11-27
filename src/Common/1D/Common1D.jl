# function ElectricField(sim::Simulation{1},m::Integer=1)
    # val = zeros(ComplexF64,3length(sim),m)
    # return ElectricField(sim,val)
# end

function ElectricField(sim::Simulation{1},val::AbstractVecOrMat)
    bl1, bl2 = sim.boundary.bls
    start = sim.x[1]
    stop = sim.x[end]
    typeof(bl1)<:AbstractComplexBL{1} ? start = start + (bl1.depth+8sim.dx) : nothing
    typeof(bl2)<:AbstractComplexBL{2} ? stop  = stop  - (bl2.depth+8sim.dx) : nothing
    start_ind = ceil(Int,Lattices.latticeindex(sim.lattice,start))
    stop_ind = floor(Int,Lattices.latticeindex(sim.lattice,stop))
    return ElectricField(sim.x,val,sim.lattice[start_ind],sim.lattice[stop_ind],[start_ind+1],[stop_ind+1])
end

lattice_primitives(constants::NTuple{1}) = (Point(1.0),), SMatrix{1,1}(1.0)

Base.getindex(lat::Lattice{1},i::CartesianIndex) = getindex(lat,i[1])
Base.getindex(lat::Lattice{1},i) = map(i->lat[i],i)
function Base.getindex(lat::Lattice{1},i::Real)
    if lat.type==:Cartesian
        return lat.origin + i*lat.dx*lat.e1
    else
        r = lat.r0 + i*lat.dr
        return lat.origin + r
    end
end

# latticeindex(lat::Lattice{1},p::Point{1}) = (p-lat.origin).x/lat.dx
latticeindex(lat::Lattice{1},p::Point{1}) = (p[1]-lat.origin[1])/lat.dx

function Base.show(io::IO, lat::Lattice{1})
    if lat.type==:Cartesian
        if !get(io, :sub, false)
            print(io, "Cartesian Lattice: \n",
                "\tdx = ", fmt("3.2f",lat.dx), "\n",
                "\torigin: ", fmt("3.2f",lat.x0))
        elseif !get(io, :sub1, false)
            print(io,
                "\n\tdx = ", fmt("3.2f",lat.dx), "\n",
                "\torigin: ", fmt("3.2f",lat.x0))
        else
            print(io,
                "\n\t\tdx = ", fmt("3.2f",lat.dx),"\n",
                "\t\torigin: ", fmt("3.2f",lat.x0))
        end
    elseif lat.type==:Polar
        println(io,"not yet implemented")
    else
        println(io, "unknown lattice type")
    end
end

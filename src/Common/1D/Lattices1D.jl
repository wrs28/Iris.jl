function Base.propertynames(::Lattice{1},private=false)
    if private
        return fieldnames(Lattice)
    else
        return (:dx, :x0)
    end
end

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
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    if lat.type==:Cartesian
        printstyled(io,"Cartesian Lattice ",color=PRINTED_COLOR_DARK)
        print(io,"(dx=")
        printstyled(io, fmt("3.4f",lat.dx),color=PRINTED_COLOR_NUMBER)
        print(io, ", origin=")
        printstyled(io, fmt("3.4f",lat.x0),color=PRINTED_COLOR_NUMBER)
        print(io,")")
    elseif lat.type==:Polar
        printstyled(io,"not yet implemented",color=PRINTED_COLOR_WARN)
    else
        printstyled(io, "unknown lattice type",color=PRINTED_COLOR_WARN)
    end
end

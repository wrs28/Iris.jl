_lattice_primitives(::TYPE,angles,constants::NTuple{1}) where TYPE = (Point{TYPE}(constants[1]),)

Base.getindex(lat::Lattice{1},i::CartesianIndex) = getindex(lat,i[1])
Base.getindex(lat::Lattice{1},i) = map(i->lat[i],i)

Base.getindex(lat::Lattice{1,Bravais} ,i::Real)     = lat.origin + i*lat.a*lat.e1
Base.getindex(lat::Lattice{1,Cartesian} ,i::Real)   = lat.origin + i*lat.dx*lat.e1
Base.getindex(lat::Lattice{1,Polar}     ,i::Real)   = lat.origin + lat.r0 + i*lat.dr
Base.getindex(lat::Lattice{1,Spherical} ,i::Real)   = lat.origin + lat.r0 + i*lat.dr

latticeindex(lat::Lattice{1,Bravais}    ,p::Point{1}) = (p.x-lat.x0)/lat.a
latticeindex(lat::Lattice{1,Cartesian}  ,p::Point{1}) = (p.x-lat.x0)/lat.dx
latticeindex(lat::Lattice{1,Polar}      ,p::Point{1}) = (p.x-lat.x0)/lat.dr - 1/2
latticeindex(lat::Lattice{1,Spherical}  ,p::Point{1}) = (p.x-lat.x0)/lat.dr - 1/2


################################################################################
# Bravais
function Base.getproperty(lat::Lattice{1,Bravais}, sym::Symbol)
    if sym == :a
        return getfield(lat,:constants)[1]
    elseif sym == :x0
        return getfield(lat,:origin)[1]
    elseif sym == :e1
        return getfield(lat,:e)[1]
    else
        return getfield(lat,sym)
    end
end

function Base.propertynames(::Lattice{1,Bravais}, private=false)
    if private
        return fieldnames(Lattice)
    else
        return (:a, :x0)
    end
end

################################################################################
# Cartesian
function Base.getproperty(lat::Lattice{1,Cartesian}, sym::Symbol)
    if sym == :dx
        return getfield(lat,:constants)[1]
    elseif sym == :x0
        return getfield(lat,:origin)[1]
    elseif sym == :e1
        return getfield(lat,:e)[1]
    else
        return getfield(lat,sym)
    end
end

function Base.propertynames(::Lattice{1,Cartesian},private=false)
    if private
        return fieldnames(Lattice)
    else
        return (:dx, :origin, :x0)
    end
end

################################################################################
# Polar
function Base.getproperty(lat::Lattice{1,Polar}, sym::Symbol)
    if sym == :dr
        return getfield(lat,:constants)[1]
    elseif sym == :r0
        return getfield(lat,:constants)[1]/2
    else
        return getfield(lat,sym)
    end
end

function Base.propertynames(::Lattice{1,Polar},private=false)
    if private
        return fieldnames(Lattice)
    else
        return (:dr, :r0)
    end
end

################################################################################
# Spherical
function Base.getproperty(lat::Lattice{1,Spherical}, sym::Symbol)
    if sym == :dr
        return getfield(lat,:constants)[1]
    elseif sym == :r0
        return getfield(lat,:constants)[1]/2
    else
        return getfield(lat,sym)
    end
end

function Base.propertynames(::Lattice{1,Spherical},private=false)
    if private
        return fieldnames(Lattice)
    else
        return (:dr, :r0)
    end
end


################################################################################
# PRETTY PRINTING
import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_WARN
import ..PRINTED_COLOR_DARK

function Base.show(io::IO, lat::Lattice{1,TYPE}) where TYPE
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    if TYPE <: Bravais
        printstyled(io,"1D Bravais Lattice ",color=PRINTED_COLOR_DARK)
        print(io,"(a=")
        printstyled(io, fmt("3.4f",lat.a),color=PRINTED_COLOR_NUMBER)
    elseif TYPE <: Cartesian
        printstyled(io,"1D Cartesian Lattice ",color=PRINTED_COLOR_DARK)
        print(io,"(dx=")
        printstyled(io, fmt("3.4f",lat.dx),color=PRINTED_COLOR_NUMBER)
    elseif TYPE <: Polar
        printstyled(io,"1D Polar Lattice ",color=PRINTED_COLOR_DARK)
        print(io,"(dr=")
        printstyled(io, fmt("3.4f",lat.dr),color=PRINTED_COLOR_NUMBER)
    elseif TYPE <: Spherical
        printstyled(io,"1D Spherical Lattice ",color=PRINTED_COLOR_DARK)
        print(io,"(dr=")
        printstyled(io, fmt("3.4f",lat.dr),color=PRINTED_COLOR_NUMBER)
    end
    print(io, ", origin=")
    printstyled(io, fmt("3.4f",lat.origin[1]),color=PRINTED_COLOR_NUMBER)
    print(io,")")
end

################################################################################
# PLOTTING
@recipe function f(lat::Lattice{1},n1=-10,n2=10)
    lat[n1:n2]
end

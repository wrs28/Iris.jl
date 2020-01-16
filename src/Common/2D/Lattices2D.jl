function _lattice_primitives(::Bravais, φ, constants::NTuple{2})
    sinφ, cosφ = sincos(φ)
    return constants[1]*Point(cosφ,sinφ), constants[2]*Point(-sinφ,cosφ)
end
function _lattice_primitives(::Cartesian, φ, constants::NTuple{2})
    sinφ, cosφ = sincos(φ[1])
    return constants[1]*Point{Cartesian}(cosφ,sinφ), constants[2]*Point{Cartesian}(-sinφ,cosφ)
end
function _lattice_primitives(::Polar, φ, constants::NTuple{2})
    return Point{Polar}(constants[1],φ[1]), Point{Polar}(constants[2],φ[1]+π/2)
end
function _lattice_primitives(::Spherical, φ, constants::NTuple{2})
    return Point{Spherical}(constants[1],φ[1]), Point{Spherical}(constants[2],φ[1]+π/2)
end


Base.getindex(lat::Lattice{2},i::CartesianIndex) = getindex(lat,i[1],i[2])
Base.getindex(lat::Lattice{2},ij::Tuple) = getindex(lat,ij[1],ij[2])
Base.getindex(lat::Lattice{2},i,j) = map((i,j)->lat[i,j],i,j)

Base.getindex(lat::Lattice{2,Bravais}   ,i::Real, j::Real)  = lat.origin + i*lat.a*lat.e1 + j*lat.b*lat.e2
Base.getindex(lat::Lattice{2,Cartesian} ,i::Real, j::Real)  = lat.origin + i*lat.dx*lat.e1 + j*lat.dy*lat.e2
Base.getindex(lat::Lattice{2,Polar}     ,i::Real, j::Real)  = lat.origin + Point{Polar}(lat.r0 + i*lat.dr, j*lat.dϕ+lat.ϕ0)
Base.getindex(lat::Lattice{2,Spherical}     ,i::Real, j::Real)  = lat.origin + Point{Spherical}(lat.r0 + i*lat.dr, j*lat.dθ+lat.θ0)

function latticeindex(lat::Lattice{2,C}, p::Point{2}) where C<:Union{Bravais,Cartesian}
    p = p-lat.origin
    x = lat.e1[1]*lat.e2[1] + lat.e1[2]*lat.e2[2]
    u1 = (lat.e1 - x*lat.e2)/(1-x^2)
    u2 = (lat.e2 - x*lat.e1)/(1-x^2)
    if C<:Bravais
        return (u1[1]*p.x+u1[2]*p.y)/lat.a, (u2[1]*p.x+u2[2]*p.y)/lat.b
    else
        return (u1[1]*p.x+u1[2]*p.y)/lat.dx, (u2[1]*p.x+u2[2]*p.y)/lat.dy
    end
end

function latticeindex(lat::Lattice{2,Polar}, p::Point{2})
    p = p-lat.origin
    return p.r/lat.dr - 1/2, (p.ϕ-lat.ϕ0)/lat.dϕ
end

function latticeindex(lat::Lattice{2,Spherical}, p::Point{2})
    p = Cartesian(p-lat.origin)
    return p.r/lat.dr - 1/2, (π/2 - p.ϕ-lat.θ0)/lat.dθ
end

################################################################################
# Bravais
function Base.getproperty(lat::Lattice{2,Bravais}, sym::Symbol)
    if sym == :a
        return getfield(lat,:constants)[1]
    elseif sym == :b
        return getfield(lat,:constants)[2]
    elseif sym == :e1
        return getfield(lat,:e)[1]
    elseif sym == :e2
        return getfield(lat,:e)[2]
    else
        return getfield(lat,sym)
    end
end

function Base.propertynames(::Lattice{2,Bravais}, private=false)
    if private
        return fieldnames(Lattice)
    else
        return (:a, :b, :origin, :e1, :e2)
    end
end

################################################################################
# Cartesian
function Base.getproperty(lat::Lattice{2,Cartesian}, sym::Symbol)
    if sym == :dx
        return getfield(lat,:constants)[1]
    elseif sym == :dy
        return getfield(lat,:constants)[2]
    elseif sym == :x0
        return getfield(lat,:origin)[1]
    elseif sym == :y0
        return getfield(lat,:origin)[2]
    elseif sym == :e1
        return getfield(lat,:e)[1]
    elseif sym == :e2
        return getfield(lat,:e)[2]
    elseif Base.sym_in(sym,(:ϕ,:φ,:phi,:𝜑,:angle))
        return getfield(lat,:e)[1].ϕ
    else
        return getfield(lat,sym)
    end
end

function Base.propertynames(::Lattice{2,Cartesian},private=false)
    if private
        return fieldnames(Lattice)
    else
        return (:dx, :dy, :e1, :e2, :origin, :ϕ)
    end
end

################################################################################
# Polar
function Base.getproperty(lat::Lattice{2,Polar}, sym::Symbol)
    if sym == :dr
        return getfield(lat,:constants)[1]
    elseif Base.sym_in(sym,(:dϕ,:dφ,:dphi))
        return getfield(lat,:constants)[2]
    elseif sym == :r0
        return getfield(lat,:constants)[1]/2
    elseif Base.sym_in(sym,(:ϕ0,:φ0,:phi0,:ϕ,:φ,:𝜑,:𝜑0,:phi,:angle))
        return getfield(lat,:e)[1].ϕ
    else
        return getfield(lat,sym)
    end
end

function Base.propertynames(::Lattice{2,Polar}, private=false)
    if private
        return fieldnames(Lattice)
    else
        return (:dr, :dϕ, :r0, :ϕ0, :origin)
    end
end


################################################################################
# Spherical
function Base.getproperty(lat::Lattice{2,Spherical}, sym::Symbol)
    if sym == :dr
        return getfield(lat,:constants)[1]
    elseif Base.sym_in(sym,(:dθ,:dϑ,:d𝜃,:dtheta))
        return getfield(lat,:constants)[2]
    elseif sym == :r0
        return getfield(lat,:constants)[1]/2
    elseif Base.sym_in(sym,(:θ0,:ϑ0,:theta0,:θ,:𝜃,:𝜃0,:ϑ,:angle))
        return acsc(getfield(lat,:e)[1][2]/getfield(lat,:e)[1][1])
    else
        return getfield(lat,sym)
    end
end

function Base.propertynames(::Lattice{2,Spherical}, private=false)
    if private
        return fieldnames(Lattice)
    else
        return (:dr, :dθ, :r0, :θ0, :origin)
    end
end



################################################################################
# PRETTY PRINTING
import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_WARN
import ..PRINTED_COLOR_DARK

function Base.show(io::IO, lat::Lattice{2,TYPE}) where TYPE
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    if TYPE <: Bravais
        printstyled(io,"2D Bravais Lattice ",color=PRINTED_COLOR_DARK)
        print(io,"(a=")
        printstyled(io, fmt("3.4f",lat.a),color=PRINTED_COLOR_NUMBER)
        print(io,", b=")
        printstyled(io, fmt("3.4f",lat.b),color=PRINTED_COLOR_NUMBER)
    elseif TYPE <: Cartesian
        printstyled(io,"2D Cartesian Lattice ",color=PRINTED_COLOR_DARK)
        print(io,"(dx=")
        printstyled(io, fmt("3.4f",lat.dx),color=PRINTED_COLOR_NUMBER)
        print(io,", dy=")
        printstyled(io, fmt("3.4f",lat.dy),color=PRINTED_COLOR_NUMBER)
    elseif TYPE <: Polar
        printstyled(io,"2D Polar Lattice ",color=PRINTED_COLOR_DARK)
        print(io,"(dr=")
        printstyled(io, fmt("3.4f",lat.dr),color=PRINTED_COLOR_NUMBER)
        print(io,", dϕ=")
        printstyled(io, fmt("3.4f",lat.dϕ),color=PRINTED_COLOR_NUMBER)
    elseif TYPE <: Spherical
        printstyled(io,"2D Spherical Lattice (w/ azimuthal symmetry) ",color=PRINTED_COLOR_DARK)
        print(io,"(dr=")
        printstyled(io, fmt("3.4f",lat.dr),color=PRINTED_COLOR_NUMBER)
        print(io,", dθ=")
        printstyled(io, fmt("3.4f",lat.dθ),color=PRINTED_COLOR_NUMBER)
    end
    print(io, ", origin=")
    printstyled(io, lat.origin, color=PRINTED_COLOR_NUMBER)
    if !iszero(lat.angle)
        print(io,", angle=")
        printstyled(io, fmt("3.2f",180*lat.angle[1]/π),"°",color=PRINTED_COLOR_NUMBER)
    end
    print(io,")")
end

################################################################################
# PLOTTING

import ..SHAPE_COLOR
import ..SHAPE_FILL_ALPHA

@recipe function f(lat::Lattice{2,Cartesian},n1=-2,n2=2)
    (lat, n1:n2, n1:n2)
end
@recipe function f(lat::Lattice{2,Cartesian},nx::AbstractVector,ny::AbstractVector)
    aspectratio --> 1
    lat[[i for i ∈ nx, j ∈ ny],[j for i ∈ nx, j ∈ ny]][:]
end
@recipe function f(lat::Lattice{2,Polar},nr=5,nϕ=floor(Int,2π/lat.dϕ))
    (lat,0:nr,0:nϕ)
end
@recipe function f(lat::Lattice{2,Polar},nr::AbstractVector,nϕ::AbstractVector)
    aspectratio --> 1
    lat[[i for i ∈ nr, j ∈ nϕ],[j for i ∈ nr, j ∈ nϕ]][:]
end
@recipe function f(lat::Lattice{2,Spherical},nr=5,nθ=floor(Int,2π/lat.dθ))
    (lat,0:nr,0:nθ)
end
@recipe function f(lat::Lattice{2,Spherical},nr::AbstractVector,nθ::AbstractVector)
    aspectratio --> 1
    lat[[i for i ∈ nr, j ∈ nθ],[j for i ∈ nr, j ∈ nθ]][:]
end
@recipe function f(lat::Lattice{2,Bravais},n1=-2,n2=2)
    (lat, n1:n2, n1:n2)
end
@recipe function f(lat::Lattice{2,Bravais},na::AbstractVector,nb::AbstractVector)
    @series begin
        aspectratio --> 1
        lat[[i for i ∈ na, j ∈ nb],[j for i ∈ na, j ∈ nb]][:]
    end
    @series begin
        alpha --> 0
        seriestype --> :path
        fillcolor --> SHAPE_COLOR
        fillrange --> mean(lat[na,nb]).y
        fillalpha --> SHAPE_FILL_ALPHA
        aspect_ratio --> 1
        legend --> false
        na0 = round(mean(na))
        nb0 = round(mean(nb))
        [lat[na0,nb0],lat[na0+1,nb0],lat[na0+1,nb0+1],lat[na0,nb0+1],lat[na0,nb0]]
    end
end



# R = @SMatrix [  cosφ    -sinφ;
                # sinφ    cosφ]
# Rᵀ = @SMatrix [  cosφ    sinφ;
#                 -sinφ    cosφ]
# return (Point(R[:,1].data),Point(R[:,2].data)), R, Rᵀ

"""
    module Lattice

for defining discrete Cartesian and Polar lattice grids used in Iris.
"""
module Lattices

# using ...Defaults
using ..Points
using Formatting
# using LinearAlgebra
# using RecipesBase
using StaticArrays

export Lattice
export latticeindex

struct Lattice{N,M,P}
    constants::NTuple{N,Float64}
    e::NTuple{N,Point{N}}
    type::Symbol
    origin::Point{N}
    angles::SVector{M,Float64}
    r0::Float64 # radius of first lattice point for Polar/Spherical
    R::SMatrix{N,N,Float64,P}

    function Lattice(constants::NTuple{N,Float64},
                type::Symbol,
                origin::Point{N},
                angles::SVector{M,Float64},
                r0::Float64) where {N,M}

        1≤N≤3 || throw(ErrorException("lattice only defined for dimensions ≤ 3, d=$N"))
        M==(N-1)N÷2 || throw(ErrorException("incorrect number of angles for dimension, should be $((N-1)N÷2), but is $M"))
        e,R = lattice_primitives(angles...,constants)
        return new{N,M,N*N}(constants,e,type,origin,angles,r0,R)
    end
end
Lattice(dx::Tuple; kwargs...) = Lattice(dx...; kwargs...)
Lattice(dx...; kwargs...) = Lattice(float.(dx); kwargs...)
function Lattice(constants::NTuple{N,Float64};
            type::Symbol = :Cartesian,
            origin = Point(ntuple(i->0.0,N)),
            angles = SVector{(N-1)N÷2,Float64}(ntuple(i->0.0,(N-1)N÷2)),
            r0::Real = constants[1]/2
            ) where N

    type = _lattice_type(type)
    origin = _lattice_origin(origin,constants)
    n = length(origin)
    (n-1)n÷2 == length(angles) || throw(ErrorException("incorrect number of angles for dimension, should be $((n-1)n÷2), but is $(length(angles))"))
    angles = _lattice_angles(angles,constants)
    return Lattice(constants,type,Point(origin),angles,r0)
end

function _lattice_type(type::Symbol)
    if type ∈ (:r,:R,:rectangular, :Rectangular, :rect, :Rect, :rectilinear, :Rectilinear, :cart, :Cart, :cartesian, :Cartesian)
        return :Cartesian
    elseif type ∈ (:p,:P,:polar, :Polar, :cyl, :Cyl, :cylindrical, :Cylindrical,:pol,:Pol)
        return :Polar
    elseif type ∈ (:s,:S,:spherical, :Spherical, :sph, :Sph)
        return :Spherical
    else
        throw("unrecognized coordinate type $type")
    end
end
_lattice_origin(origin,constants::NTuple{N}) where N = Point(origin)
_lattice_angles(angles,constants::NTuple{N}) where N = SVector{(N-1)N÷2,Float64}(angles)

include("1D/Lattices.jl")
# include("2D/Lattices.jl")
# include("3D/Lattices.jl")

function Base.getproperty(lat::Lattice,sym::Symbol)
    if sym == :dx
        return getfield(lat,:constants)[1]
    elseif sym == :dy
        return getfield(lat,:constants)[2]
    elseif sym == :dz
        return getfield(lat,:constants)[3]
    elseif sym == :dr
        return getfield(lat,:constants)[1]
    elseif sym == :x0
        return getfield(lat,:origin)[1]
    elseif sym == :y0
        return getfield(lat,:origin)[2]
    elseif sym == :z0
        return getfield(lat,:origin)[3]
    elseif sym == :e1
        return getfield(lat,:e)[1]
    elseif sym == :e2
        return getfield(lat,:e)[2]
    elseif sym == :e3
        return getfield(lat,:e)[3]
    elseif Base.sym_in(sym,(:dϕ,:dφ,:dphi))
        return getfield(lat,:constants)[2]
    elseif Base.sym_in(sym,(:dθ,:dϑ,:dtheta))
        return getfield(lat,:constants)[3]
    elseif Base.sym_in(sym,(:ϕ,:φ,:phi))
        return getfield(lat,:angles)[1]
    elseif Base.sym_in(sym,(:θ,:ϑ,:theta))
        return getfield(lat,:angles)[1]
    else
        return getfield(lat,sym)
    end
end

latticeindex(lat::Lattice,x...) = latticeindex(lat,Point(x...))

end # module

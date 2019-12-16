"""
for defining discrete Cartesian and Polar lattice grids used in Iris.
"""
module Lattices

export Lattice
export Bravais
export latticeindex

files = (
    "1D/Lattices1D.jl",
    "2D/Lattices2D.jl",
    # "3D/Lattices3D.jl"
)

using ..Points
using Formatting
using LinearAlgebra
using RecipesBase
using StaticArrays
using Statistics

import ..Points.Cartesian
import ..Points.Polar
import ..Points.Spherical
AbstractLatticeType = Points.AbstractCoordinateType
struct Bravais <: AbstractLatticeType end

struct Lattice{N,TYPE,CE,CO}
    constants::NTuple{N,Float64} # lattice spacings
    e::NTuple{N,Point{N,CE}} # unit vectors
    origin::Point{N,CO}

    function Lattice(::TYPE,
                primitives::NTuple{N,Point{N,CE}},
                origin::Point{N,CO},
                ) where {N,CE,CO} where TYPE<:AbstractLatticeType

        1≤N≤3 || throw(ErrorException("lattice only defined for dimensions ≤ 3, d=$N"))
        if TYPE<:Polar
            2≤N≤3 || throw(ErrorException("lattice type $TYPE not consistent with dimension $N, which must be ≥2"))
        end
        if TYPE<:Spherical
            3≤N≤3 || throw(ErrorException("lattice type $TYPE not consistent with dimension $N, which must be 3"))
        end
        constants = _lattice_constants(primitives)
        vectors = _lattice_vectors(primitives)
        return new{N,TYPE,CE,CO}(constants, vectors, origin)
    end
end

_lattice_constants(primitives::Tuple) = map(p->norm(p),primitives)
_lattice_vectors(primitives::Tuple) = map(p->p/norm(p),primitives)

"""
    Bravais(primitives...; [origin]) -> bravaislattice
"""
Bravais(args...; kwargs...) = Lattice(Bravais(), args...; kwargs...)

"""
    Cartesian(dx,...; kwargs) -> clat

    Cartesian((dx,...); kwargs...) -> clat
"""
Cartesian(args...; kwargs...) = Lattice(Cartesian(), args...; kwargs...)

"""
    Polar(dx,...; [origin]) -> polarlat
"""
Polar(args...; kwargs...) = Lattice(Polar(),args...; kwargs...)

"""
    Spherical(dx,...; [origin]) -> sphlat
"""
Spherical(dx...; kwargs...) = Lattice(Spherical(), dx...; kwargs...)

Lattice(t,dx::Tuple; kwargs...) = Lattice(t,dx...; kwargs...)
Lattice(t,dx::Number,dy...; kwargs...) = Lattice(t,(float(dx),float.(dy)...); kwargs...)

"""
    Lattice(::Bravais, primitives...; [origin]) -> bravaislattice
"""
Lattice(::Bravais, primitive::Point,primitives...; kwargs... ) = Lattice(Bravais(),(primitive,primitives...);kwargs...)
function Lattice(::Bravais,
            primitives::NTuple{N,Point};
            origin = Point(ntuple(i->0.0,N)),
            ) where N
    constants = _lattice_constants(primitives)
    return Lattice(Bravais(), Cartesian.(primitives), _lattice_origin(origin,constants))
end
"""
    Lattice(::Union{Cartesian,Polar,Spherical}, dx...; [origin]) -> lattice
"""
function Lattice(type::Union{Cartesian,Polar,Spherical},
            constants::NTuple{N,Float64};
            origin = Point(ntuple(i->0.0,N)),
            angles = zeros(Float64,(N-1)N÷2),#SVector{M,Float64}(ntuple(i->0.0,(N-1)N÷2)),
            ) where N

    M = length(angles)
    M==(N-1)N÷2 || throw(ErrorException("incorrect number of angles for dimension, should be $((N-1)N÷2), but is $M"))
    return Lattice(type, _lattice_primitives(type,angles,constants), _lattice_origin(origin,constants))
end


_lattice_origin(origin::Point{N},constants::NTuple{N}) where N = origin
_lattice_origin(origin,constants::NTuple{N}) where N = Point(origin...)


foreach(include,files)

"""
    latticeindex(lattice, point) -> ind::Float64
"""
latticeindex(lat::Lattice,x::Real,y...) = latticeindex(lat,Point(x,y...))

end # module

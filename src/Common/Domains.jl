"""
for defining Lattice, Nondispersive, and Dispersive Domains
"""
module Domains

files = (
    "1D/Domains1D.jl",
    "2D/Domains2D.jl",
    # "3D/Domains3D.jl"
)

export AbstractDomain
export LatticeDomain
export NondispersiveDomain
export DispersiveDomain
export Symmetric
export Unsymmetric

using ..Boundaries
using ..DielectricFunctions
using ..PumpFunctions
using ..Dispersions
using ..Lattices
using ..Points
using ..Shapes
using LinearAlgebra
using RecipesBase

"""
	Unsymmetric
"""
struct Unsymmetric end
"""
	Symmetric
"""
LinearAlgebra.Symmetric

"""
    AbstractDomain{N}
"""
abstract type AbstractDomain{N} end

################################################################################
# LATTICE DOMAIN

struct LatticeDomain{N,CLASS,TC,TBND,TLAT,TP} <: AbstractDomain{N}
    boundary::TBND
    lattice::TLAT
    n::ComplexF64
    ε::ComplexF64
    type::Symbol
    name::Symbol
    x::Vector{Point{N,TP}}
    indices::Vector{CartesianIndex{N}}
    interior::BitArray{1}
    bulk::BitArray{1}
    surface::BitArray{1}
    corner::BitArray{1}
    nnm::NTuple{N,Vector{Int}}
    nnp::NTuple{N,Vector{Int}}

    function LatticeDomain{CLASS}(
        boundary::TBND,
        lattice::TLAT,
        n::Number,
        type::Symbol,
        name::Symbol,
        x::Vector{Point{N,TP}},
        indices::Vector{CartesianIndex{N}},
        interior::BitArray{1},
        bulk::BitArray{1},
        surface::BitArray{1},
        corner::BitArray{1},
        nnm::NTuple{N,Vector{Int}},
        nnp::NTuple{N,Vector{Int}}
        ) where {TBND<:Boundary{N},TP,TLAT<:Lattice{N,TC},CLASS} where {N,TC}

        new{N,CLASS,TC,TBND,TLAT,TP}(boundary,lattice,n,n^2,type,name,x,indices,interior,bulk,surface,corner,nnm,nnp)
    end
end

function Base.getproperty(ldom::LatticeDomain, sym::Symbol)
    if sym == :shape
        return getfield(getfield(ldom,:boundary),:shape)
    elseif Base.sym_in(sym, propertynames(getfield(ldom,:lattice)))
        return getproperty(getfield(ldom,:lattice),sym)
    else
        return getfield(ldom,sym)
    end
end

################################################################################
# NONDISPERSIVE DOMAIN

struct NondispersiveDomain{N,TSH,TDF} <: AbstractDomain{N}
    shape::TSH
    dielectric::TDF
    type::Symbol
    name::Symbol

    function NondispersiveDomain(
                shape::TSH,
                dielectric::TDF,
                type::Symbol = :generic,
                name::Symbol = :anonymous
                ) where {TSH<:AbstractShape{N}, TDF<:AbstractDielectricFunction} where N

        new{N,TSH,TDF}(shape,dielectric,type,name)
    end
end

NondispersiveDomain(shape::AbstractShape, n::Number=1, args...) =
    NondispersiveDomain(shape,DielectricFunctions.PiecewiseConstant(n), args...)

NondispersiveDomain(shape::AbstractShape, n1::Number, n2::Number, args...) =
    NondispersiveDomain(shape,DielectricFunctions.PiecewiseConstant(n1,n2), args...)

NondispersiveDomain(boundary::Boundary, args...) = NondispersiveDomain(boundary.shape, args...)

################################################################################
# DISPERSIVE DOMAIN

struct DispersiveDomain{N,TSH,TCHI,TPF} <: AbstractDomain{N}
    shape::TSH
    χ::TCHI
    pump::TPF
    type::Symbol
    name::Symbol

    function DispersiveDomain(
                shape::TSH,
                χ::TCHI,
                pump::TPF,
                type::Symbol = :generic,
                name::Symbol = :anonymous
                ) where {TSH<:AbstractShape{N}, TPF<:AbstractPumpFunction, TCHI<:AbstractDispersion} where N

        new{N,TSH,TCHI,TPF}(shape,χ,pump,type,name)
    end
end

DispersiveDomain(shape::AbstractShape, χ::AbstractDispersion, F::Number=1, args...) =
    DispersiveDomain(shape, χ, PumpFunctions.PiecewiseConstant(F), args...)

DispersiveDomain(boundary::Boundary, args...) = DispersiveDomain(boundary.shape, args...)


# load dimensional files
foreach(include,files)


################################################################################
# Pretty Printing

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK

function Base.show(io::IO,dom::LatticeDomain{N,CLASS}) where {N,CLASS}
    print(io,"$(N)D $CLASS ")
    printstyled(io,"LatticeDomain ",color=PRINTED_COLOR_DARK)
    println(io,"(",dom.type,"): ",dom.name)
    print(io,"\tNumber of sites: ")
    printstyled(io,length(dom.x),"\n",color=PRINTED_COLOR_NUMBER)
    print(io,"\tNumber of bulk points: ")
    printstyled(io,sum(dom.bulk),"\n",color=PRINTED_COLOR_NUMBER)
    print(io,"\tNumber of surface points: ")
    printstyled(io,sum(dom.surface),"\n",color=PRINTED_COLOR_NUMBER)
    print(io,"\tNumber of corner points: ")
    printstyled(io,sum(dom.corner),"\n",color=PRINTED_COLOR_NUMBER)
    print(io,"\t\t====================")
    println(io)
    println(IOContext(io,:tabbed2=>true),dom.lattice)
    println(io)
    print(IOContext(io,:tabbed2=>true),dom.boundary)
end

function Base.show(io::IO,ndom::NondispersiveDomain)
    printstyled(io,"NondispersiveDomain ",color=PRINTED_COLOR_DARK)
    println(io,"(",ndom.type,"): ",ndom.name)
    println(io,"\t\t",ndom.shape)
    show(IOContext(io,:tabbed2=>true), ndom.dielectric)
end

function Base.show(io::IO,ddom::DispersiveDomain)
    printstyled(io,"DispersiveDomain ",color=PRINTED_COLOR_DARK)
    println(io,"(",ddom.type,"): ",ddom.name)
    println(io,"\t\t",ddom.shape)
    show(IOContext(io,:tabbed2=>true),ddom.pump)
    println(io)
    println(io,"\t\t",ddom.χ)
end

end # module

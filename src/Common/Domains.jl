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

using ..Boundaries
using ..DielectricFunctions
using ..PumpFunctions
using ..Dispersions
using ..Lattices
using ..Points
using ..Shapes
using RecipesBase

import ..Symmetric, ..Unsymmetric
import LinearAlgebra.norm

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
        x::Vector{Point{N,TP}},
        indices::Vector{CartesianIndex{N}},
        interior::BitArray{1},
        bulk::BitArray{1},
        surface::BitArray{1},
        corner::BitArray{1},
        nnm::NTuple{N,Vector{Int}},
        nnp::NTuple{N,Vector{Int}};
        type::Symbol = :generic,
        name::Symbol = :anonymous
        ) where {TBND<:Boundary{N},TP,TLAT<:Lattice{N,TC},CLASS} where {N,TC}

        new{N,CLASS,TC,TBND,TLAT,TP}(boundary,lattice,n,n^2,type,name,x,indices,interior,bulk,surface,corner,nnm,nnp)
    end
end

LatticeDomain(lattice::Lattice, boundary::Boundary, n::Number=1; kwargs...) =
    LatticeDomain(boundary, lattice, n; kwargs...)

LatticeDomain(boundary::Boundary, n::Number, lattice::Lattice; kwargs...) =
    LatticeDomain(boundary, lattice, n; kwargs...)

LatticeDomain(lattice::Lattice, n::Number, boundary::Boundary; kwargs...) =
    LatticeDomain(boundary, lattice, n; kwargs...)

LatticeDomain(n::Number, boundary::Boundary, lattice::Lattice; kwargs...) =
    LatticeDomain(boundary, lattice, n; kwargs...)

LatticeDomain(n::Number, lattice::Lattice, boundary::Boundary; kwargs...) =
    LatticeDomain(boundary, lattice, n; kwargs...)

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


Base.propertynames(ndom::NondispersiveDomain, private=false) = (fieldnames(NondispersiveDomain)..., propertynames(ndom.dielectric,private)...)

function Base.getproperty(ndom::NondispersiveDomain, sym::Symbol)
    dielectric = getfield(ndom,:dielectric)
    if Base.sym_in(sym,propertynames(dielectric,true))
        return getproperty(dielectric,sym)
    else
        return getfield(ndom,sym)
    end
end

function Base.setproperty!(ndom::NondispersiveDomain, sym::Symbol, val)
    if Base.sym_in(sym,propertynames(ndom.dielectric,true))
        return setproperty!(ndom.dielectric,sym,val)
    else
        return setfield!(ndom,sym,val)
    end
end

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
                pump::TPF;
                type::Symbol = :generic,
                name::Symbol = :anonymous
                ) where {TSH<:AbstractShape{N}, TPF<:AbstractPumpFunction, TCHI<:AbstractDispersion} where N

        new{N,TSH,TCHI,TPF}(shape,χ,pump,type,name)
    end
end

DispersiveDomain(shape::AbstractShape, χ::AbstractDispersion, F::Number=1; kwargs...) =
    DispersiveDomain(shape, χ, PumpFunctions.PiecewiseConstant(F); kwargs...)

DispersiveDomain(boundary::Boundary, args...; kwargs...) = DispersiveDomain(boundary.shape, args...; kwargs...)


Base.propertynames(ddom::DispersiveDomain, private=false) = (fieldnames(DispersiveDomain)..., propertynames(ddom.pump,private)...)

function Base.getproperty(ddom::DispersiveDomain, sym::Symbol)
    pump = getfield(ddom,:pump)
    if Base.sym_in(sym,propertynames(pump,true))
        return getproperty(pump,sym)
    else
        return getfield(ddom,sym)
    end
end

function Base.setproperty!(ddom::DispersiveDomain, sym::Symbol, val)
    if Base.sym_in(sym,propertynames(ddom.pump,true))
        return setproperty!(ddom.pump,sym,val)
    else
        return setfield!(ddom,sym,val)
    end
end

################################################################################
# load dimensional files (mostly lattice domain)
foreach(include,files)


################################################################################
# Pretty Printing

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_VARIABLE

function Base.show(io::IO,dom::LatticeDomain{N,CLASS}) where {N,CLASS}
    print(io,"$(N)D ")
    CLASS<:Symmetric ? print(io,"Symmetric ") : nothing
    CLASS<:Unsymmetric ? print(io,"Unsymmetric ") : nothing
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
    print(io,"\tBackground Index: ")
    printstyled(io,dom.n;color=PRINTED_COLOR_NUMBER)
    println(io)
    print(io,"\t\t====================")
    println(io)
    println(IOContext(io,:tabbed2=>true),dom.lattice)
    println(io)
    print(IOContext(io,:tabbed2=>true),dom.boundary)
end

function Base.show(io::IO,ndom::NondispersiveDomain)
    printstyled(io,"NondispersiveDomain ",color=PRINTED_COLOR_DARK)
    println(io,"(",ndom.type,"): ",ndom.name)
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    println(io,"\t\t",ndom.shape)
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    show(IOContext(io,:tabbed2=>true), ndom.dielectric)
end

function Base.show(io::IO,ddom::DispersiveDomain)
    printstyled(io,"DispersiveDomain ",color=PRINTED_COLOR_DARK)
    println(io,"(",ddom.type,"): ",ddom.name)
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    println(io,"\t\t",ddom.shape)
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    show(IOContext(io,:tabbed2=>true),ddom.pump)
    println(io)
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    println(io,"\t\t",ddom.χ)
end

end # module

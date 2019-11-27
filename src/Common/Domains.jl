"""
    module Domains

for constructing nondispersive domains to be fed to Simulation
"""
module Domains

files = (
    "1D/Domains1D.jl",
    # "2D/Domains2D.jl",
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

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK

abstract type AbstractDomain{N} end

"""
    struct LatticeDomains

    NondispersiveDomain(type,boundary,lattice,[dielectric,pump];align=false,[name]) -> domain
    NondispersiveDomain(type,lattice,boundary,[dielectric,pump];align=false,[name]) -> domain

combines `boundary` with `lattice` (in either order) to generate a list of sites.
`T` just labels the kinds of domain (e.g. Cavity, or Waveguide)
each site is labeled by being in the interior, or containing the interior (`domain.surface`),
or being a corner.

----------------------
    (::NondispersiveDomain)(args...; align=false) -> dom

construct a new domain with modified parameters in `args`
This is for "updating" non-geometric fields of the immutable `domain`.
For example, to change the polarity of the boundary layers, do
    `new_dom = old_dom(conj(old_dom.boundary))`

For geometric parameters it recomputes the whole domain, so this command is no
more efficient than explicitly constructing a new domain. It might save some typing, however.
"""
struct LatticeDomain{N,TBND,TLAT} <: AbstractDomain{N}
    boundary::TBND
    lattice::TLAT
    n::ComplexF64
    ε::ComplexF64
    type::Symbol
    name::Symbol
    x::Vector{Point{N}}
    indices::Vector{CartesianIndex{N}}
    interior::BitArray{1}
    bulk::BitArray{1}
    surface::BitArray{1}
    corner::BitArray{1}
    nnm::NTuple{N,Vector{Int}}
    nnp::NTuple{N,Vector{Int}}

    function LatticeDomain(boundary::TBND, lattice::TLAT, n::Number, type::Symbol,
        name::Symbol, x::Vector{Point{N}}, indices::Vector{CartesianIndex{N}},
        interior::BitArray{1}, bulk::BitArray{1}, surface::BitArray{1}, corner::BitArray{1},
        nnm::NTuple{N,Vector{Int}}, nnp::NTuple{N,Vector{Int}}
        ) where {TBND<:Boundary{N},TLAT<:Lattice{N}} where N

        new{N,TBND,TLAT}(boundary,lattice,n,n^2,type,name,x,indices,interior,bulk,surface,corner,nnm,nnp)
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

function Base.propertynames(::LatticeDomain, private=false)
    if private
        return fieldnames(LatticeDomain)
    else
        return (:boundary, :lattice, :shape, :n, :ε, :type, :name, :x)
    end
end

"""
    struct NondispersiveDomain

    NondispersiveDomain(type,boundary,lattice,[dielectric,pump];align=false,[name]) -> domain
    NondispersiveDomain(type,lattice,boundary,[dielectric,pump];align=false,[name]) -> domain

combines `boundary` with `lattice` (in either order) to generate a list of sites.
`T` just labels the kinds of domain (e.g. Cavity, or Waveguide)
each site is labeled by being in the interior, or containing the interior (`domain.surface`),
or being a corner.

----------------------
    (::NondispersiveDomain)(args...; align=false) -> dom

construct a new domain with modified parameters in `args`
This is for "updating" non-geometric fields of the immutable `domain`.
For example, to change the polarity of the boundary layers, do
    `new_dom = old_dom(conj(old_dom.boundary))`

For geometric parameters it recomputes the whole domain, so this command is no
more efficient than explicitly constructing a new domain. It might save some typing, however.
"""
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


"""
    struct DispersiveDomain

    DispersiveDomain(type,boundary,lattice,[dielectric,pump];align=false,[name]) -> domain
    DispersiveDomain(type,lattice,boundary,[dielectric,pump];align=false,[name]) -> domain

combines `boundary` with `lattice` (in either order) to generate a list of sites.
`T` just labels the kinds of domain (e.g. Cavity, or Waveguide)
each site is labeled by being in the interior, or containing the interior (`domain.surface`),
or being a corner.

----------------------
    (::DispersiveDomain)(args...;align=false) -> dom

construct a new domain with modified parameters in `args`
This is for "updating" non-geometric fields of the immutable `domain`.
For example, to change the polarity of the boundary layers, do
    `new_dom = old_dom(conj(old_dom.boundary))`

For geometric parameters it recomputes the whole domain, so this command is no
more efficient than explicitly constructing a new domain. It might save some typing, however.
"""
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

function Base.show(io::IO,dom::LatticeDomain)
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




# @recipe function f(d::Domain)
#     aspect_ratio --> 1
#     legend --> false
#     @series begin
#         seriestype --> :scatter
#         markersize --> MARKERSIZE_SCALE/sqrt(length(d.x))
#         markerstrokealpha --> 0
#         shape --> MARKERSHAPE
#         (d.x[d.bulk],d.y[d.bulk])
#     end
#     @series begin
#         seriestype --> :scatter
#         markersize --> MARKERSIZE_SCALE/sqrt(length(d.x))
#         markerstrokealpha --> 0
#         shape --> MARKERSHAPE
#         (d.x[d.surface],d.y[d.surface])
#     end
#     @series begin
#         seriestype --> :scatter
#         markersize --> MARKERSIZE_SCALE/sqrt(length(d.x))
#         markerstrokealpha --> 0
#         shape --> MARKERSHAPE
#         (d.x[d.corner],d.y[d.corner])
#     end
#     @series d.boundary
# end
#








# function (dom::Domain{N})(dielectric::DielectricFunction) where N
#     return new{N,typeof(bnd.boundary),typeof(bnd.lattice),typeof(bnd.dielectric),typeof(bnd.pump)}(dom.type,
#     dom.name,dom.boundary,dom.lattice,dielectric,dom.pump,
#     dom.x,dom.indices,dom.ε,dom.F,dom.interior,dom.bulk,
#     dom.surface,dom.corner,dom.nnm,dom.nnp)
# end
# function (dom::Domain{N})(pump::PumpFunction) where N
#     return newnew{N,typeof(bnd.boundary),typeof(bnd.lattice),typeof(bnd.dielectric),typeof(bnd.pump)}(dom.type,
#     dom.name,dom.boundary,dom.lattice,dom.dielectric,pump,
#     dom.x,dom.indices,dom.ε,dom.F,dom.interior,dom.bulk,
#     dom.surface,dom.corner,dom.nnm,dom.nnp)
# end
# function (dom::Domain{N})(boundary::Boundary) where N
#     return newnew{N,typeof(bnd.boundary),typeof(bnd.lattice),typeof(bnd.dielectric),typeof(bnd.pump)}(dom.type,
#     dom.name,boundary,dom.lattice,dom.dielectric,dom.pump,
#     dom.x,dom.indices,dom.ε,dom.F,dom.interior,dom.bulk,
#     dom.surface,dom.corner,dom.nnm,dom.nnp)
# end






# # for modifying non-geometric parameters of an already existing domain
# function (dom::Domain)(args...)
#     for a ∈ args
#         dom = dom(a)
#     end
#     return dom
# end

"""
    module Domains

for constructing domains to be fed to Simulation
"""
module Domains

files = (
    "1D/Domains1D.jl",
    # "2D/Domains2D.jl",
    # "3D/Domains3D.jl"
)

export Domain

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK
using ..Boundaries
using ..DielectricFunctions
using ..Dispersions
using ..Lattices
using ..Points
using ..Shapes


"""
    struct Domain

    Domain(type,boundary,lattice,[dielectric,pump];align=false,[name]) -> domain
    Domain(type,lattice,boundary,[dielectric,pump];align=false,[name]) -> domain

combines `boundary` with `lattice` (in either order) to generate a list of sites.
`T` just labels the kinds of domain (e.g. Cavity, or Waveguide)
each site is labeled by being in the interior, or containing the interior (`domain.surface`),
or being a corner.

----------------------
    (::Domain)(args...;align=false) -> dom

construct a new domain with modified parameters in `args`
This is for "updating" non-geometric fields of the immutable `domain`.
For example, to change the polarity of the boundary layers, do
    `new_dom = old_dom(conj(old_dom.boundary))`

For geometric parameters it recomputes the whole domain, so this command is no
more efficient than explicitly constructing a new domain. It might save some typing, however.
"""
struct Domain{N,TBND,TLAT,TDF,TPF,TCHI}
    type::Symbol
    name::Symbol
    boundary::TBND
    lattice::TLAT
    dielectric::TDF
    pump::TPF

    x::Array{Point{N},1}
    indices::Array{CartesianIndex{N},1}
    ε::Array{ComplexF64,1}
    F::Array{Float64,1}
    χ::TCHI

    interior::BitArray{1}
    bulk::BitArray{1}
    surface::BitArray{1}
    corner::BitArray{1}

    nnm::NTuple{N,Array{Int,1}}
    nnp::NTuple{N,Array{Int,1}}

    function Domain(
        type::Symbol,
        name::Symbol,
        boundary::TBND,
        lattice::TLAT,
        dielectric::TDF,
        pump::TPF,
        x::Array{Point{N},1},
        indices::Array{CartesianIndex{N},1},
        ε::Array{ComplexF64,1},
        F::Array{Float64,1},
        χ::TC,
        interior::BitArray{1},
        bulk::BitArray{1},
        surface::BitArray{1},
        corner::BitArray{1},
        nnm::NTuple{N,Array{Int,1}},
        nnp::NTuple{N,Array{Int,1}}
        ) where {TBND<:Boundary,TLAT<:Lattice{N},TDF<:DielectricFunction,TPF<:PumpFunction,TC<:Tuple} where N

        new{N,TBND,TLAT,TDF,TPF,TC}(type,name,boundary,lattice,dielectric,pump,x,indices,ε,F,χ,interior,bulk,surface,corner,nnm,nnp)
    end

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
end

function Base.conj(d::Domain)
    Domain(
        d.type,
        d.name,
        conj(d.boundary),
        d.lattice,
        conj(d.dielectric),
        d.pump,
        d.x,
        d.indices,
        conj(d.ε),
        d.F,
        map(conj,d.χ),
        d.interior,
        d.bulk,
        d.surface,
        d.corner,
        d.nnm,
        d.nnp
        )
end


# # for modifying non-geometric parameters of an already existing domain
# function (dom::Domain)(args...)
#     for a ∈ args
#         dom = dom(a)
#     end
#     return dom
# end

foreach(include,files)

function Base.show(io::IO,dom::Domain)
    printstyled(io,"Domain ",color=PRINTED_COLOR_DARK)
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
    println(IOContext(io,:tabbed2=>true),dom.boundary)
    println(io)
    println(IOContext(io,:tabbed2=>true),dom.dielectric)
    println(io)
    println(IOContext(io,:tabbed2=>true),dom.pump)
    println(io)
    println(io,"\t\t",dom.χ)
end

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

end # module

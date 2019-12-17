"""
For defining boundary conditions, used in building LatticeDomains.
"""
module Boundaries

export Boundary

files = (
	"1D/Boundaries1D.jl",
	"2D/Boundaries2D.jl",
	# "3D/Boundaries3D.jl"
	)

using ..BoundaryConditions
using ..BoundaryLayers
using ..Shapes
using RecipesBase
using TupleTools

import ..BL_DEPTH
import ..DEFAULT_BC
import ..DEFAULT_BL

@eval begin
	const DEFAULT_BCS = $(DEFAULT_BC)
	const DEFAULT_BLS = $(DEFAULT_BL)
end
const DEFAULT_BL_DEPTH = BL_DEPTH

"""
	Boundary{N}

`N`-dimensional boundary
Fields: `shape`, `bcs`, `bls`
"""
struct Boundary{N,TS,TBC,TBL}
    shape::TS
    bcs::TBC
    bls::TBL

    function Boundary(
				shape::AbstractShape{NDIMS,NSIDES},
		        bcs::NTuple{NSIDES,AbstractBC},
		        bls::NTuple{NSIDES,AbstractBL}
		        ) where {NDIMS,NSIDES}

		length(bcs)==nsides(shape) || throw("incorrect number of boundary conditions $(length(bcs)), should match number of sides $(nsides(sh))")
		length(bls)==nsides(shape) || throw("incorrect number of boundary layers $(length(bls)), should match number of sides $(nsides(sh))")

		bcs_sorted = map(i->bcs[i],TupleTools.sortperm(map(b->getside(b),bcs)))
		bls_sorted = map(i->bls[i],TupleTools.sortperm(map(b->getside(b),bls)))
		bls_sorted = map(x->_bls_by_shape(x,shape),bls_sorted)
	    return new{NDIMS,typeof(shape),typeof(bcs_sorted),typeof(bls_sorted)}(shape,bcs_sorted,bls_sorted)
    end
end

_bls_by_shape(bls::AbstractBL,shape::AbstractShape) = bls

foreach(include,files)

"""
	Boundary(shape, boundary_conditions..., boundary_layers...; depth=$DEFAULT_BL_DEPTH) -> bnd

`boundary_conditions` can be given as a tuple, individual conditions (such as
`DirichletBC{2}()` for Dirichlet on side 2), or a type to be applied to all sides,
such as [`NeumannBC`](@ref). Any side left otherwise unspecified is assumed to be `$DEFAULT_BCS`.
The boundary condition types are [`DirichletBC`](@ref), [`NeumannBC`](@ref), [`MatchedBC`](@ref), [`FloquetBC`](@ref),
`noBC`.

`boundary_layers` is specified similarly to `boundary_conditions`. The types are
[`PML`](@ref), [`cPML`](@ref), [`noBL`](@ref) (defaults to [`$DEFAULT_BLS`](@ref)). The optional argument `depth` sets the depth
of any unspecified layer (defaults to `$DEFAULT_BL_DEPTH`).

No default `shape` is assumed.

The arguments can be entered in any order. For example, `Boundary(PML{2}(.2),shape,NeumannBC)`,
`Boundary((NeumannBC{1}(),DirichletBC{2}()),shape,cPML)` and `Boundary(DirichletBC,shape)` are
all valid.
"""
function Boundary(
			shape::AbstractShape{NDIMS,NSIDES},
			bcs::NTuple{MBC,AbstractBC},
			bls::NTuple{MBL,AbstractBL}
			) where {NDIMS,NSIDES,MBC,MBL}

	MBC ≤ NSIDES || throw("more boundary conditions $MBC than sides $NSIDES")
	MBL ≤ NSIDES || throw("more boundary layers $MBC than sides $NSIDES")

	all_sides = ntuple(identity,NSIDES)

	sides_bcs = map(b->getside(b),bcs)
	if isempty(sides_bcs)
		missing_sides_bcs = all_sides
	else
		missing_sides_bcs = TupleTools.deleteat(all_sides,sides_bcs)
	end
	new_bcs = map(s->DEFAULT_BCS{s}(),missing_sides_bcs)

	sides_bls = map(b->getside(b),bls)
	if isempty(sides_bls)
		missing_sides_bls = all_sides
	else
		missing_sides_bls = TupleTools.deleteat(all_sides,sides_bls)
	end
	new_bls = map(s->DEFAULT_BLS{s}(DEFAULT_BL_DEPTH),missing_sides_bls)

	return Boundary(shape,(bcs...,new_bcs...), (bls...,new_bls...))
end

"""
	ndims(::Boundary) -> n
"""
Base.ndims(bnd::Boundary) = ndims(bnd.shape)

"""
	nsides(::Boundary) -> n
"""
Shapes.nsides(bnd::Boundary) = nsides(bnd.shape)

function Boundary(args...;depth::Real=DEFAULT_BL_DEPTH)
	sh = _getshape(args...)
	bcs = _getbcs(sh,args...)
	bls = _getbls(sh,args...)
	for bl ∈ bls
		bl ∈ args ? nothing : bl.depth=depth
	end
	return Boundary(sh,bcs,bls)
end

_getshape(sh::AbstractShape,args...) = sh
_getshape(arg,args...) = _getshape(args...)
_getshape() = throw("no shape given")

_getbcs(sh::AbstractShape,bc::AbstractBC,args...) = (bc,_getbcs(sh::AbstractShape,args...)...)
_getbcs(sh::AbstractShape,::Type{BC},args...) where BC<:AbstractBC = ntuple(i->BC{i}(),nsides(sh))
_getbcs(sh::AbstractShape,arg,args...) = _getbcs(sh::AbstractShape,args...)
_getbcs(sh::AbstractShape) = ()

_getbls(sh::AbstractShape,bl::AbstractBL,args...) = (bl,_getbls(sh,args...)...)
_getbls(sh::AbstractShape,::Type{BL},args...) where BL<:AbstractBL = ntuple(i->BL{i}(DEFAULT_BL_DEPTH),nsides(sh))
_getbls(sh::AbstractShape,arg,args...) = _getbls(sh,args...)
_getbls(sh::AbstractShape) = ()

################################################################################
# PRETTY PRINTING

import ..PRINTED_COLOR_DARK

function Base.show(io::IO,bnd::Boundary)
	get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
	printstyled(io,"Boundary\n",color=PRINTED_COLOR_DARK)
	get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
	print(io,"\t",bnd.shape)
	for i ∈ 1:nsides(bnd)
		println(io)
		get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
		print(io,"\tSIDE ",i,":")
		print(IOContext(io,:compact=>true),"\t", bnd.bcs[i],"  &  ",bnd.bls[i])
	end
end

end # module

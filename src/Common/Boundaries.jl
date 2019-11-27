"""
	module Boundaries
"""
module Boundaries

export Boundary

files = (
	"1D/Boundaries1D.jl",
	# "2D/Boundaries2D.jl",
	# "3D/Boundaries3D.jl"
	)

using ..BoundaryConditions
using ..BoundaryLayers
using ..Shapes
using TupleTools

import ..PRINTED_COLOR_DARK
import ..BL_DEPTH
import ..DEFAULT_BC
import ..DEFAULT_BL

@eval begin
	const DEFAULT_BCS = $(DEFAULT_BC)
	const DEFAULT_BLS = $(DEFAULT_BL)
end
const DEFAULT_BL_DEPTH = BL_DEPTH

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

		bcs_sorted = map(i->bcs[i],TupleTools.sortperm(map(b->get_side(b),bcs)))
		bls_sorted = map(i->bls[i],TupleTools.sortperm(map(b->get_side(b),bls)))
		bls_sorted = map(x->_bls_by_shape(x,shape),bls_sorted)
	    return new{NDIMS,typeof(shape),typeof(bcs_sorted),typeof(bls_sorted)}(shape,bcs_sorted,bls_sorted)
    end
end

foreach(include,files)

"""
	Boundary(shape, boundary_conditions..., boundary_layers...; depth=$DEFAULT_BL_DEPTH) -> bnd

`boundary_conditions` can be given as a tuple, individual conditions (such as
`DirichletBC{2}()` for Dirichlet on side 2), or a type to be applied to all sides,
such as `NeumannBC`. Any side left otherwise unspecified is assumed to be `noBC`.
The boundary condition types are `DirichletBC`, `NeumannBC`, `MatchedBC`, `FloquetBC`,
`noBC` (the default).

`boundary_layers` is specified similarly to `boundary_conditions`. The types are
`PML`, `cPML`, `noBL` (the default). The optional argument `depth` sets the depth
of any unspecified layer.

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

	sides_bcs = map(b->get_side(b),bcs)
	if isempty(sides_bcs)
		missing_sides_bcs = all_sides
	else
		missing_sides_bcs = TupleTools.deleteat(all_sides,sides_bcs)
	end
	new_bcs = map(s->DEFAULT_BCS{s}(),missing_sides_bcs)

	sides_bls = map(b->get_side(b),bls)
	if isempty(sides_bls)
		missing_sides_bls = all_sides
	else
		missing_sides_bls = TupleTools.deleteat(all_sides,sides_bls)
	end
	new_bls = map(s->DEFAULT_BLS{s}(DEFAULT_BL_DEPTH),missing_sides_bls)

	return Boundary(shape,(bcs...,new_bcs...), (bls...,new_bls...))
end

Base.ndims(bnd::Boundary) = ndims(bnd.shape)
Shapes.nsides(bnd::Boundary) = nsides(bnd.shape)

function Boundary(args...;depth::Real=DEFAULT_BL_DEPTH)
	sh = _get_shape(args...)
	bcs = _get_bcs(sh,args...)
	bls = _get_bls(sh,args...)
	for bl ∈ bls
		bl ∈ args ? nothing : bl.depth=depth
	end
	return Boundary(sh,bcs,bls)
end

_get_shape(sh::AbstractShape,args...) = sh
_get_shape(arg,args...) = _get_shape(args...)
_get_shape() = throw("no shape given")

_get_bcs(sh::AbstractShape,bc::AbstractBC,args...) = (bc,_get_bcs(sh::AbstractShape,args...)...)
_get_bcs(sh::AbstractShape,::Type{BC},args...) where BC<:AbstractBC = ntuple(i->BC{i}(),nsides(sh))
_get_bcs(sh::AbstractShape,arg,args...) = _get_bcs(sh::AbstractShape,args...)
_get_bcs(sh::AbstractShape) = ()

_get_bls(sh::AbstractShape,bl::AbstractBL,args...) = (bl,_get_bls(sh,args...)...)
_get_bls(sh::AbstractShape,::Type{BL},args...) where BL<:AbstractBL = ntuple(i->BL{i}(DEFAULT_BL_DEPTH),nsides(sh))
_get_bls(sh::AbstractShape,arg,args...) = _get_bls(sh,args...)
_get_bls(sh::AbstractShape) = ()


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


# """
	# conj(::Boundary) -> bnd
#
# Make a new `Boundary` object with conjugated boundary layers.
# """
# Base.conj(bnd::Boundary) = Boundary(bnd.shape,conj.(bnd.bcs),conj.(bnd.bls))

_bls_by_shape(bls,shape::AbstractShape) = bls

end # module

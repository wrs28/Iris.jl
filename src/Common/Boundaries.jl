"""
	module Boundaries
"""
module Boundaries

files = (
	"1D/Boundaries1D.jl",
	# "2D/Boundaries2D.jl",
	# "3D/Boundaries3D.jl"
	)

export Boundary

import ..PRINTED_COLOR_DARK

# using ...Defaults
using ..BoundaryConditions
using ..BoundaryLayers
using ..Shapes
using TupleTools

const DEFAULT_BC = noBC
const DEFAULT_BL = noBL

"""
	struct Boundary

	Boundary(shape::Shape,bcs,bls) -> bnd

`bcs` is either tuple of boundary conditions, or it is a single boundary condition which applies to all sides.
`bls` is either tuple of boundary layers, or it is a single boundary layer which applies to all sides.

if fewer boundary conditions or layers are supplied than the number of sides of
`shape`, it assumes `noBC` and `noBL`
"""
struct Boundary{TS,TBC,TBL}
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
	    return new{typeof(shape),typeof(bcs_sorted),typeof(bls_sorted)}(shape,bcs_sorted,bls_sorted)
    end
end


function Boundary(
	shape::AbstractShape{NDIMS,NSIDES},
	bcs::NTuple{MBC,AbstractBC},
	bls::NTuple{MBL,AbstractBC}
	) where {NDIMS,NSIDES,MBC,MBL}

	MBC ≤ NDIMS || throw("more boundary conditions $MBC than sides $NSIDES")
	MBL ≤ NDIMS || throw("more boundary layers $MBC than sides $NSIDES")

	all_sides = ntuple(identity,NSIDES)

	sides_bcs = map(b->get_side(b),bcs)
	missing_sides_bcs = _get_missing_sides(all_sides,sides_bcs)
	new_bcs = ntuple(i->DEFAULT_BC{missing_sides_bcs[i]}(),NSIDES-MBC)

	sides_bls = map(b->get_side(b),bls)
	missing_sides_bls = _get_missing_sides(all_sides,sides_bls)
	new_bls = ntuple(i->DEFAULT_BL{missing_sides_bls[i]}(),NSIDES-MBL)

	return Boundary(shape,(bcs...,new_bcs...), (bls...,new_bls...))
end

function _get_missing_sides(all_sides::NTuple,sides::NTuple)
	if all_sides[1] ∈ sides
		return _get_missing_sides(Base.tail(all_sides),sides)
	else
		return (all_sides[1],_get_missing_sides(Base.tail(all_sides),sides)...)
	end
end
_get_missing_sides(all_sides::Tuple{},sides::NTuple) = ()

Base.ndims(bnd::Boundary) = ndims(bnd.shape)
Shapes.nsides(bnd::Boundary) = nsides(bnd.shape)

function Boundary(sh::AbstractShape,::Type{BC}=DEFAULT_BC,::Type{BL}=DEFAULT_BL; depth::Real=.1) where {BC<:AbstractBC,BL<:AbstractBL}
	bcs = ntuple(i->BC{i}(),nsides(sh))
	bls = ntuple(i->BL{i}(1),nsides(sh))
	for bl ∈ bls bl.depth = depth end
	return Boundary(sh,bcs,bls)
end
Boundary(sh::AbstractShape,bl::Type{BL},bc::Type{BC}=DEFAULT_BC; depth::Real=.1) where {BC<:AbstractBC,BL<:AbstractBL} =
	Boundary(sh,bc,bl;depth=depth)
Boundary(bc::Type{BC},sh::AbstractShape,bl::Type{BL}=DEFAULT_BL; depth::Real=.1) where {BC<:AbstractBC,BL<:AbstractBL} =
	Boundary(sh,bc,bl;depth=depth)
Boundary(bl::Type{BL},sh::AbstractShape,bc::Type{BC}=DEFAULT_BC; depth::Real=.1) where {BC<:AbstractBC,BL<:AbstractBL} =
	Boundary(sh,bc,bl;depth=depth)
Boundary(bc::Type{BC},bl::Type{BL},sh::AbstractShape; depth::Real=.1) where {BC<:AbstractBC,BL<:AbstractBL} =
	Boundary(sh,bc,bl;depth=depth)
Boundary(bl::Type{BL},bc::Type{BC},sh::AbstractShape; depth::Real=.1) where {BC<:AbstractBC,BL<:AbstractBL} =
	Boundary(sh,bc,bl;depth=depth)


function Boundary(args...)
	sh = _get_shape(args...)
	bcs = _get_bcs(args...)
	bls = _get_bls(args...)
	return Boundary(sh,bcs,bls)
end
_get_shape(sh::AbstractShape,args...) = sh
_get_shape(args...) = _get_shape(Base.tail(args)...)
_get_shape() = throw("no shape given")

_get_bcs(bc::AbstractBC,args...) = (bc,_get_bcs(args...)...)
_get_bcs(args...) = _get_bcs(Base.tail(args)...)
_get_bcs() = ()

_get_bls(bc::AbstractBL,args...) = (bc,_get_bls(args...)...)
_get_bls(args...) = _get_bls(Base.tail(args)...)
_get_bls() = ()


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


"""
	conj(::Boundary) -> bnd

Make a new `Boundary` object with conjugated boundary layers.
"""
Base.conj(bnd::Boundary) = Boundary(bnd.shape,conj.(bnd.bcs),conj.(bnd.bls))

_bls_by_shape(bls,shape::AbstractShape) = bls

foreach(include,files)

end # module

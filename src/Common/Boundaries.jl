"""
	module Boundaries
"""
module Boundaries

export Boundary

using ...Defaults
using ..BoundaryConditions
using ..BoundaryLayers
using ..Shapes

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

    function Boundary(shape::AbstractShape{NDIMS,NSIDES},
        bcs::NTuple{NSIDES,AbstractBC},
        bls::NTuple{NSIDES,AbstractBL}
        ) where {NDIMS,NSIDES}

		bls = map(x->_bls_by_shape(x,shape),bls)
	    return new{typeof(shape),typeof(bcs),typeof(bls)}(shape,bcs,bls)
    end
end

Base.ndims(bnd::Boundary) = ndims(bnd.shape)
Shapes.nsides(bnd::Boundary) = nsides(bnd.shape)


# function Boundary(shape::AbstractShape{NS},bc::AbstractBC=noBC(),bls::AbstractBL=noBL()) where N
	# bcs = fill(bcs,N)
	# bls = fill(bls,N)
	# return Boundary(shape,(bcs...,),(bls...,))
# end
# function Boundary(shape::AbstractShape{N},bcs::AbstractBC,bls) where N
# 	bcs = fill(bcs,N)
# 	return Boundary(shape,(bcs...,),bls)
# end
# function Boundary(shape::AbstractShape{N},bcs,bls::AbstractBL=noBL()) where N
# 	bls = fill(bls,N)
# 	return Boundary(shape,bcs,(bls...,))
# end


function Base.show(io::IO,bnd::Boundary)
	N = nsides(bnd)
	print(io,"Boundary: \n")
	print(io,"\tshape: $(bnd.shape)")
	for i ∈ 1:N
		println(io)
		print(io,"\t", bnd.bcs[i],"  +  ",bnd.bls[i])
	end
end


"""
	conj(::Boundary) -> bnd

Make a new `Boundary` object with conjugated boundary layers.
"""
Base.conj(bnd::Boundary) = Boundary(bnd.shape,conj.(bnd.bcs),conj.(bnd.bls))

_bls_by_shape(bls,shape::AbstractShape) = bls

include("1D/Boundaries.jl")
# include("2D/Boundaries.jl")
# include("3D/Boundaries.jl")

end # module

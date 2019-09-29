"""
    module Shapes

Definitions of various shapes.

All shapes have the signature: `shape(parameters::Tuple, x0, y0, θ=0)`, and some
have keyword argument `reference` which determines which point `(x0,y0)` is referring to.

The length of `parameters` depends on the shape. For example, `Circle` has just one
parameter, the radius, while `Ellipse` has two, the semi-major and semi-minor axes.
"""
module Shapes

export AbstractShape
export nsides

abstract type AbstractShape{NDIMS,NSIDES} end
Base.ndims(::AbstractShape{NDIMS,NSIDES}) where {NDIMS,NSIDES} = NDIMS
nsides(::AbstractShape{NDIMS,NSIDES}) where {NDIMS,NSIDES} = NSIDES

include("1D/Shapes.jl")
# include("2D/Shapes.jl")
# include("3D/Shapes.jl")

end #module

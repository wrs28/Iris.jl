"""
Definitions of various shapes in 1D, 2D, 3D.

All shapes have the signature: `shape(params::Tuple, x0, y0, θ=0)`, and some
have keyword argument `reference` which determines which point `(x0,y0)` is referring to.

The length of `parameters` depends on the shape. For example, `Circle` has just one
parameter, the radius, while `Ellipse` has two, the semi-major and semi-minor axes.
"""
module Shapes

files = (
    "1D/Shapes1D.jl",
    "2D/Shapes2D.jl",
    # "3D/Shapes3D.jl",
    )

export AbstractShape
export nsides

using ..Points

"""
    AbstractShape{NDIMS,NSIDES}
"""
abstract type AbstractShape{NDIMS,NSIDES} end

"""
    ndims(shape) = d
"""
Base.ndims(::AbstractShape{NDIMS,NSIDES}) where {NDIMS,NSIDES} = NDIMS

"""
    nsides(shape) = n
"""
nsides(::AbstractShape{NDIMS,NSIDES}) where {NDIMS,NSIDES} = NSIDES

foreach(include,files)

end #module

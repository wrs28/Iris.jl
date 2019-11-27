#TODO: make ε a tensor
"""
    module PumpFunctions

for defining pump functions either as piecewise constant
or with user-supplied functions

To define a user-supplied PumpFunction, define a `struct` which is a subtype
of `AbstractPumpFunction`, and implement `(::Userdefined)(p::Point)` whose
output is a real scalar, the pump.
"""
module PumpFunctions

export AbstractPumpFunction

using ..Points

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK

abstract type AbstractPumpFunction <: Function end

"""
    PiecewiseConstant([F=1]) -> pc
"""
mutable struct PiecewiseConstant <: AbstractPumpFunction F::Float64 end
PiecewiseConstant() = PiecewiseConstant(1)

"""
    PiecewiseConstant(point) -> F
"""
(pc::PiecewiseConstant)(args...) = pc.F

# Base.conj(pc::PiecewiseConstant) = pc


# Pretty Printing
function Base.show(io::IO,d::PiecewiseConstant)
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    printstyled(io, "PiecewiseConstant",color=PRINTED_COLOR_DARK)
    print(io, " (")
    printstyled(io,"F = ",d.F,color=PRINTED_COLOR_NUMBER)
    print(io,")")
end
Base.show(io::IO,::MIME"text/plain",d::PiecewiseConstant) = show(io,d)

end # module

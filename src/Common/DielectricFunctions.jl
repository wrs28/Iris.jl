#TODO: make ε a tensor
"""
    module DielectricFunctions

for defining dielectric functions either as piecewise constant
or with user-supplied functions

To define a user-supplied DielectricFunction, define a `struct` which is a subtype
of `AbstractDielectricFunction`, and implement `(::Userdefined)(p::Point)` whose
output is a complex scalar, the dielectric.
"""
module DielectricFunctions

export AbstractDielectricFunction

using ..Points

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK

abstract type AbstractDielectricFunction <: Function end

"""
    PiecewiseConstant([n1=1,n=0]) -> pc

    PiecewiseConstant(n::Complex) -> pc
"""
mutable struct PiecewiseConstant <: AbstractDielectricFunction
    n₁::Float64
    n₂::Float64
    ε::ComplexF64
    PiecewiseConstant(n1::Real=1,n2::Real=0) = new(n1,n2,complex(n1,n2)^2)
end
PiecewiseConstant(n::Complex) = PiecewiseConstant(real(n),imag(n))

"""
    PiecewiseConstant(point) -> ε
"""
(pc::PiecewiseConstant)(p::Point) = pc.ε

function Base.getproperty(pc::PiecewiseConstant, sym::Symbol)
    if sym == :n1
        return getfield(pc,:n₁)
    elseif sym == :n2
        return getfield(pc,:n₂)
    else
        return getfield(pc,sym)
    end
end

function Base.setproperty!(pc::PiecewiseConstant, sym::Symbol, x)
    if sym == :n1
        setfield!(pc,:n₁,x)
    elseif sym == :n2
        setfield!(pc,:n₂,x)
    elseif sym == :ε
        n1,n2 = reim(sqrt(x))
        setfield!(pc,:n₁,n1)
        setfield!(pc,:n₂,n2)
    else
        setfield!(pc,sym,x)
    end
    setfield!(pc,:ε,complex(getfield(pc,:n₁),getfield(pc,:n₂))^2)
    return nothing
end

# Base.conj(pc::DF) = PiecewiseConstant(pc.n1,-de.n2)

# Pretty Printing
function Base.show(io::IO,d::PiecewiseConstant)
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    printstyled(io, "PiecewiseConstant",color=PRINTED_COLOR_DARK)
    print(io, " (")
    printstyled(io,"n = ",d.n₁+1im*d.n₂,color=PRINTED_COLOR_NUMBER)
    print(io,")")
end
Base.show(io::IO,::MIME"text/plain",d::PiecewiseConstant) = show(io,d)

end # module

"""
    module Dispersions
"""
module Dispersions

export AbstractDispersion
export NoDispersion
export TwoLevelSystem
export susceptability

abstract type AbstractDispersion end
struct NoDispersion <: AbstractDispersion end
susceptability(::NoDispersion, args...) = complex(0.0)

include("Dispersions/TwoLevelSystems.jl")
using .TwoLevelSystems

include("Dispersions/Kerrs.jl")
using .Kerrs

end # module

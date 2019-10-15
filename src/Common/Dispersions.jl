"""
    module Dispersions
"""
module Dispersions

export AbstractDispersion
export TwoLevelSystem
export susceptability
export jacobian_lasing

abstract type AbstractDispersion end

susceptability(χs::Tuple,args...) = map(x->susceptability(x,args...),χs)

include("Dispersions/TwoLevelSystems.jl")
using .TwoLevelSystems

include("Dispersions/Kerrs.jl")
using .Kerrs

end # module

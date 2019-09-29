"""
    module Dispersions
"""
module Dispersions

export AbstractDispersion
export TwoLevelSystem

abstract type AbstractDispersion end

include("Dispersions/TwoLevelSystem.jl")
include("Dispersions/Kerr.jl")

end # module

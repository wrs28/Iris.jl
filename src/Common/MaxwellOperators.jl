"""
    module MaxwellOperators
"""
module MaxwellOperators

export maxwell
export maxwell_lep
export maxwell_nep

# # using ...Defaults
# # import ..Shapes
using ..Simulations
# # using Base.Iterators
using LinearAlgebra
using NonlinearEigenproblems
using SparseArrays

include("1D/MaxwellOperators.jl")
# include("2D/MaxwellOperators.jl")
# include("3D/MaxwellOperators.jl")

end #module

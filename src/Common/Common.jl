"""
Essential tools used throughout Iris.

Submodules:
  * [`Boundaries`](@ref)
  * [`BoundaryLayers`](@ref)
  * [`BoundaryConditions`](@ref)
  * [`Curlcurls`](@ref)
  * [`DielectricFunctions`](@ref)
  * [`Dispersions`](@ref)
  * [`Domains`](@ref)
  * [`ElectricFields`](@ref)
  * [`Lattices`](@ref)
  * [`LU_Factorizations`](@ref)
  * [`Plotting`](@ref)
  * [`Points`](@ref)
  * [`PumpFunctions`](@ref)
  * [`SelfEnergies`](@ref)
  * [`Shapes`](@ref)
  * [`Simulations`](@ref)
"""
module Common

dimensional_files = (
    "1D/Common1D.jl",
)

# `Defaults`: contains all default parameters in one location for uniformity and easy alteration.
include("../Defaults.jl")

# Base.conj(::Tuple{}) = ()

include("Points.jl")
using .Points
export Point
export Cartesian
export Polar
export Spherical

include("VectorFields.jl")
using .VectorFields
export ScalarField
export ElectricField
export VectorField

include("Shapes.jl")
using .Shapes
export AbstractShape
export Interval
export Circle
export Square
export Rectangle

include("BoundaryLayers.jl")
using .BoundaryLayers
export PML
export cPML
export noBL

include("BoundaryConditions.jl")
using .BoundaryConditions
export noBC
export DirichletBC
export NeumannBC
export FloquetBC
export MatchedBC
export BCHermiticity
export HermitianBC
export NonHermitianBC

include("Boundaries.jl")
using .Boundaries
export Boundary

include("DielectricFunctions.jl")
using .DielectricFunctions
export AbstractDielectricFunction

include("PumpFunctions.jl")
using .PumpFunctions
export AbstractPumpFunction

include("Dispersions.jl")
using .Dispersions
export TwoLevelSystem
export susceptability

include("Lattices.jl")
using .Lattices
export Lattices
export Lattice
export latticeindex

include("Domains.jl")
using .Domains
export LatticeDomain
export NondispersiveDomain
export DispersiveDomain
export Symmetric
export Unsymmetric

# include("Laplacians.jl")
# using .Laplacians
# export Laplacian

# include("Curlcurls.jl")
# using .Curlcurls
# export Curlcurl

# # include("Curls.jl")
# # using .Curls
# # export Curl

# include("SelfEnergies.jl")
# using .SelfEnergies

# include("Simulations.jl")
# using .Simulations
# export Simulation
# export smooth!
# export update!
# export update_dielectric!
# export update_pump!
# export Unsymmetric
# export Symmetric
# export Hermitian

# include("HelmholtzOperators.jl")
# using .HelmholtzOperators
# export Helmholtz

# include("MaxwellOperators.jl")
# using .MaxwellOperators
# export Maxwell

# include("LU_Factorizations.jl")
# using .LU_Factorizations
# export DEFAULT_LUPACK
# export MSolver
# export PSolver
# export USolver
# export AbstractLUPACK

# include("Plotting.jl")
# using .Plotting

# foreach(include,dimensional_files)

end # module

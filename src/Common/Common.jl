"""
Essential tools used throughout Iris.

Submodules:
  * [`Points`](@ref)
  * [`ElectricFields`](@ref)
  * [`Shapes`](@ref)
  * [`BoundaryLayers`](@ref)
  * [`BoundaryConditions`](@ref)
  * [`Boundaries`](@ref)
  * [`DielectricFunctions`](@ref)
  * [`PumpFunctions`](@ref)
  * [`Dispersions`](@ref)
  * [`Lattices`](@ref)
  * [`Domains`](@ref)
  * [`Curlcurls`](@ref)
  * [`SelfEnergies`](@ref)
  * [`Simulations`](@ref)
  * [`LU_Factorizations`](@ref)
  * [`Plotting`](@ref)
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

include("ElectricFields.jl")
using .ElectricFields
export ElectricField

include("Shapes.jl")
using .Shapes
export AbstractShape
export Interval
export Circle
export Ellipse
export DeformedDisk
export Annulus
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

include("Domains.jl")
using .Domains
export LatticeDomain
export NondispersiveDomain
export DispersiveDomain
# export Cavity
# export Resonator
# export Void
# export Dielectric
# export Waveguide
# export Lead

include("Curlcurls.jl")
using .Curlcurls
export Curlcurl
#
# include("Curls.jl")
# using .Curls
# export Curl

include("SelfEnergies.jl")
using .SelfEnergies

include("Simulations.jl")
using .Simulations
export Simulation
export smooth!
export update!
export update_dielectric!
export update_pump!
export Unsymmetric
export Symmetric
export Hermitian

include("MaxwellOperators.jl")
using .MaxwellOperators
export Maxwell
export maxwell

include("LU_Factorizations.jl")
using .LU_Factorizations
export DEFAULT_LUPACK
export MSolver
export PSolver
export USolver
export AbstractLUPACK

include("Plotting.jl")
using .Plotting

foreach(include,dimensional_files)

end # module

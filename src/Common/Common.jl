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

"""
	Symmetric
"""
struct Symmetric end

"""
	Unsymmetric
"""
struct Unsymmetric end

# `Defaults`: contains all default parameters in one location for uniformity and easy alteration.
include("../Defaults.jl")

include("Points.jl") # ✅
using .Points
export Point
export Cartesian
export Polar
export Spherical

include("VectorFields.jl") # ✅
using .VectorFields
export ScalarField
export ElectricField
export VectorField
export update!
export component

include("Shapes.jl") # ✅
using .Shapes
export AbstractShape
export Interval
export Circle
export Square
export Rectangle

include("BoundaryLayers.jl") # ✅
using .BoundaryLayers
export PML
export cPML
export noBL

include("BoundaryConditions.jl") # ✅
using .BoundaryConditions
export noBC
export DirichletBC
export NeumannBC
export FloquetBC
export MatchedBC
export BCHermiticity
export HermitianBC
export NonHermitianBC

include("Boundaries.jl") # ✅
using .Boundaries
export Boundary

include("DielectricFunctions.jl") # ✅
using .DielectricFunctions
export AbstractDielectricFunction

include("PumpFunctions.jl") # ✅
using .PumpFunctions
export AbstractPumpFunction

include("Dispersions.jl") # ✅
using .Dispersions
export NoDispersion
export TwoLevelSystem
export susceptability

include("Lattices.jl") # ✅ Except for Bravais
using .Lattices
export Lattices
export Lattice
export latticeindex

include("Domains.jl") # ✅ Only 1D, 2D Cartesian
using .Domains
export Domains
export LatticeDomain
export NondispersiveDomain
export DispersiveDomain

include("Laplacians.jl") # ✅ Only 1D, 2D Cartesian
using .Laplacians
export Laplacian

include("Curlcurls.jl")
using .Curlcurls
export Curlcurl

# # include("Curls.jl")
# # using .Curls
# # export Curl

include("SelfEnergies.jl") # ✅ Only 1D, 2D Cartesian Symmetric Dirichlet & Open, functions need work
using .SelfEnergies

include("Simulations.jl") # ✅ Only 1D, 2D Cartesian Symmetric sans Curlcurls and no smoothing in 2D
using .Simulations
export Simulation
export smooth!
export update_dielectric!
export update_pump!

include("HelmholtzOperators.jl") # ✅ Only 1D, 2D Cartesian Symmetric
using .HelmholtzOperators
export Helmholtz

include("MaxwellOperators.jl")
using .MaxwellOperators
export Maxwell

include("LU_Factorizations.jl") # ✅ checked for UMFPACK, not for Pardiso or MUMPS
using .LU_Factorizations
export DEFAULT_LUPACK
export MSolver
export PSolver
export USolver
export AbstractLUPACK

include("Plotting.jl") # ✅ basic 2D and 1D simulation plotting
using .Plotting

end # module

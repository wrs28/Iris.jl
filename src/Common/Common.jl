"""
    module Common

Essential tools used in the rest of Iris.
"""
module Common

# `Defaults`: contains all default parameters in one location for uniformity and easy alteration.
include("../Defaults.jl")

Base.conj(::Tuple{}) = ()

include("Points.jl") # ✅
using .Points
export Point

include("ElectricFields.jl") # ✅
using .ElectricFields
export ElectricField

include("Shapes.jl") # ✅
using .Shapes
export AbstractShape
export Interval
export Circle
export Ellipse
export DeformedDisk
export Annulus
export Square
export Rectangle

include("BoundaryLayers.jl") # ✅
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

include("Boundaries.jl")
using .Boundaries
export Boundary

include("DielectricFunctions.jl")
using .DielectricFunctions
export DielectricFunction
export PumpFunction

include("Dispersions.jl")
using .Dispersions
export TwoLevelSystem
export jacobian_lasing
export susceptability

include("Lattices.jl")
using .Lattices
export Lattices
export Lattice

include("Domains.jl")
using .Domains
export Domain
export Cavity
export Resonator
export Void
export Dielectric
export Waveguide
export Lead

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
export update_dielectric!
export update_pump!
export update!
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

ElectricField(sim::Simulation{1}) = ElectricField(sim.x,1)
ElectricField(sim::Simulation{1},val) = ElectricField(sim.x,val)
ElectricField(val,sim::Simulation{1}) = ElectricField(val,sim.x)

end # module

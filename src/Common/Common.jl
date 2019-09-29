"""
    module IrosBase

Essential tools used in the rest of Iros.
"""
module Common

include("Points.jl")
using .Points
export Point

include("Shapes.jl")
using .Shapes
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

include("SelfEnergies.jl")
using .SelfEnergies

include("Simulations.jl")
using .Simulations
export Simulation
export smooth!
export update_dielectric!
export update_pump!
export Unsymmetric
export Symmetric
export Hermitian

include("MaxwellOperators.jl")
using .MaxwellOperators
export maxwell
export maxwell_lep
export maxwell_nep

end # module

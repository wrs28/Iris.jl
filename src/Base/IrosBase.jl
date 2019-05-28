"""
    module IrosBase

Essential tools used in the rest of Iros.
"""
module IrosBase

include("Shapes.jl")
using .Shapes
export Circle,
Ellipse,
DeformedDisk,
Annulus,
Square,
Rectangle,
AbstractShape

include("DielectricFunctions.jl")
using .DielectricFunctions
export DielectricFunction

include("Lattices.jl")
using .Lattices
export Lattices

include("Boundaries.jl")
using .Boundaries
export Boundary,
noBC,
DirichletBC,
NeumannBC,
FloquetBC,
MatchedBC,
PML,
cPML,
noBL

include("Tessellations.jl")
using .Tessellations

include("Domains.jl")
using .Domains
export Domain,
Cavity,
Resonator,
Void,
Dielectric,
Waveguide,
Lead

include("Simulations.jl")
using .Simulations
export Simulation,
smooth!,
update_dielectric!,
update_pump!

include("Laplacians.jl")
using .Laplacians
export Laplacian

end # module

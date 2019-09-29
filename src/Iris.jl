"""
    module Iris

A package for scalar finite difference resonance and scattering calculations, for both
    closed and open systems.
    * `Defaults`: contains all default parameters in one location for uniformity and easy alteration.
    * `IrosBase`: the essential structures used to define a simulation and to construct the necessary operators
    * `Spectral`: linear and non-linear eigenvalue problems for resonant, anti-resonant, and RSM problems,
        for closed systems, or for open systems via boundary matching or PMLs.
    * `Floquet`: for finding the band structures and dispersion curves of crystal structures
    * `Scattering`: for solving inhomogeneous scattering problems
    * `SALT`: the SALT algorithm for lasing
    * `SCPA`: saturable CPA
"""
module Iris

using Reexport

include("Defaults.jl")
using .Defaults

include("Common/Common.jl")
@reexport using .Common

include("Spectral/Spectral.jl")
@reexport using .Spectral
#
# include("Floquet/Floquet.jl")
# @reexport using .Floquet
#
# include("Scattering/Scattering.jl")
# @reexport using .Scattering
#
# include("SALT/SALT.jl")
# @reexport using .SALT
#
# include("SaturableAbsorption/SCPA.jl")
# @reexport using .SCPA

end # module

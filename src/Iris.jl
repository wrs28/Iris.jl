"""
    module Iris

A package for scalar finite difference resonance and scattering calculations, for both
    closed and open systems.

    * `Defaults`: contains all default parameters in one location for uniformity and easy alteration.
    * `Common`: the essential structures used to define a simulation and to construct the necessary operators
    * `Spectral`: linear and non-linear eigenvalue problems for resonant, anti-resonant, and RSM problems,
        for closed systems, or for open systems via boundary matching or PMLs.
    * `Floquet`: for finding the band structures and dispersion curves of crystal structures
    * `Scattering`: for solving inhomogeneous scattering problems
    * `Lasing`: the SALT algorithm for lasing
    * `SaturableCPA`: saturable CPA
    * `TimeDomain: `
"""
module Iris

using Reexport

include("Common/Common.jl")
@reexport using .Common

include("Spectral/Spectral.jl")
@reexport using .Spectral

include("Floquet/Floquet.jl")
@reexport using .Floquet

include("Scattering/Scattering.jl")
@reexport using .Scattering

include("Lasing/Lasing.jl")
@reexport using .Lasing

include("SaturableCPA/SaturableCPA.jl")
@reexport using .SaturableCPA

include("TimeDomain/TimeDomain.jl")
@reexport using .TimeDomain

end # module

"""
    module Iros

A package for scalar finite difference resonance and scattering calculations, for both
    closed and open systems.
    * `IrosDefaults`: contains all default parameters in one location for uniformity and easy alteration.
    * `IrosBase`: the essential structures used to define a simulation and to construct the necessary operators
    * `IrosSpectral`: linear and non-linear eigenvalue problems for resonant, anti-resonant, and RSM problems,
        for closed systems, or for open systems via boundary matching or PMLs.
    * `IrosFloquet`: for finding the band structures and dispersion curves of crystal structures
    * `IrosScattering`: for solving inhomogeneous scattering problems
    * `IrosSALT`: the SALT algorithm for lasing
    * `IrosSCPA`: saturable CPA
"""
module Iros

using Reexport

include("Base/IrosBase.jl")
@reexport using IrosBase

include("Spectral/IrosSpectral.jl")
@reexport using IrosSpectral

include("Floquet/IrosFloquet.jl")
@reexport using IrosFloquet

include("Scattering/IrosScattering.jl")
@reexport using IrosScattering

include("SALT/IrosSALT.jl")
@reexport using IrosSALT

include("SCPA/IrosSCPA.jl")
@reexport using IrosSCPA

end # module

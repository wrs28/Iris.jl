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

@reexport using IrosBase
@reexport using IrosSpectral
@reexport using IrosFloquet
@reexport using IrosScattering
@reexport using IrosSALT
@reexport using IrosSCPA

end # module

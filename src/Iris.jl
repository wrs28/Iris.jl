"""
    module Iris

A package for scalar finite difference resonance, scattering calculations, and
    time-domain simulations for both closed and open systems.

    * `Common`: the essential structures used to define a simulation and to construct the necessary operators
    * `Spectral`: linear and non-linear eigenvalue problems for resonant, anti-resonant, and RSM problems,
        for closed systems, or for open systems via boundary matching or PMLs.
    * `Floquet`: for finding the band structures and dispersion curves of crystal structures
    * `Scattering`: for solving inhomogeneous scattering problems, linear and nonlinear
    * `Lasing`: the SALT algorithm for lasing
    * `SaturableCPA`: the SALT algorithm for CPA
    * `TimeDomain`: time evolution of Maxwell and Maxwell-Bloch field equations
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

include("SMatrices/SMatrices.jl")
@reexport using  .SMatrices

include("Lasing/Lasing.jl") #TODO: SALTbootstrap [ ]; boundary-condition-modifying convenience wrappers; multimode Jacobian ✅
@reexport using .Lasing

include("SaturableCPA/SaturableCPA.jl") #TODO: boundary-condition-modifying convenience wrappers
@reexport using .SaturableCPA

include("TimeDomain/TimeDomain.jl")
@reexport using .TimeDomain

end # module

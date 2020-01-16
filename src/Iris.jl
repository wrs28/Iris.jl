"""
Electromagnetic resonance, scattering, and time-domain finite difference
simulations for both closed and open systems.

Divided across seven submodules:

  * [`Common`](@ref): the essential structures used to define a simulation and to construct the necessary differential operators
  * [`Spectral`](@ref): linear and non-linear eigenvalue problems for resonance, anti-resonance, and RSM problems, for closed systems, or for open systems via boundary matching or PMLs.
  * [`Floquet`](@ref): for finding the band structures and dispersion curves of crystal structures
  * [`Scattering`](@ref): for solving inhomogeneous scattering problems, linear and nonlinear
  * [`SMatrices`](@ref): for S-matrix computations (works in parallel)
  * [`Lasing`](@ref): the SALT algorithm for lasing
  * [`SaturableCPA`](@ref): the SALT algorithm for CPA
  * [`TimeDomain`](@ref): time evolution of Maxwell and Maxwell-Bloch field equations
"""
module Iris

using Reexport

include("Common/Common.jl")
@reexport using .Common
import .Common: Symmetric, Unsymmetric

include("Spectral/Spectral.jl")
@reexport using .Spectral

include("Floquet/Floquet.jl")
@reexport using .Floquet

include("Scattering/Scattering.jl")
@reexport using .Scattering

include("SMatrices/SMatrices.jl")
@reexport using  .SMatrices

include("Lasing/Lasing.jl")
@reexport using .Lasing

include("SaturableCPA/SaturableCPA.jl")
@reexport using .SaturableCPA

include("TimeDomain/TimeDomain.jl")
@reexport using .TimeDomain

end # module

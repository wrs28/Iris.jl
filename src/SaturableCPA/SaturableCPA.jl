# TODO: add simulation-modifying convenience wrappers to make bc's incoming. Without that, this module provides essentially no utility above SALT
"""
    module SaturableCPA

Convenience wrappers for solution of the nonlinear SatruableCPA equation.
Exports `HelmholtzSCPA` and `MaxwellSCPA` constructors and extends methods of
`NLsolve`, which must be separately imported.
"""
module SaturableCPA

export HelmholtzSCPA
export MaxwellSCPA

using ..Common
using ..Lasing
using NLsolve
using RecipesBase

import ..Lasing.SALTProblem
import ..Lasing.HelmholtzSALT
import ..Lasing.MaxwellSALT

################################################################################
# SCPA objects and constructors (mostly hooks into SALT constructors)

"""
    SCPAProblem{M,THS}

Saturable CPA object with `M` component fields
"""
struct SCPAProblem{NMODES,M,THS}
    SALT::THS

    SCPAProblem(SALT::SALTProblem{NMODES,M}) where {NMODES,M} = new{NMODES,M,typeof(SALT)}(SALT)
end

function Base.getproperty(mcpa::SCPAProblem, sym::Symbol)
    if sym == :SALT
        getfield(mcpa,:SALT)
    else
        getproperty(getfield(mcpa,:SALT),sym)
    end
end

Base.propertynames(::SCPAProblem,private=false) = propertynames(SALTProblem,private)


HelmholtzSCPA{NMODES} = SCPAProblem{NMODES,1}
"""
    HelmholtzSCPA(sim::Simulation, m::Integer) -> scpa

SPCA object with `m` fields for Helmholtz problem
"""
HelmholtzSCPA(sim::Simulation,m::Integer) = SCPAProblem(HelmholtzSALT{m}(sim))


MaxwellSCPA{NMODES} = SCPAProblem{NMODES,3}
"""
    MaxwellSCPA(sim::Simulation, m::Integer) -> scpa

SPCA object with `m` fields for Maxwell problem
"""
MaxwellSCPA(sim::Simulation,m::Integer) = SCPAProblem(MaxwellSALT{m}(sim))



################################################################################
# Extend NLsolve functions

fnames = (:nlsolve,:fixedpoint)
for fn ∈ fnames
    @eval NLsolve.$(fn)(ms::SCPAProblem; kwargs...) = $(fn)(ms, ms.ωs, ms.ψs; kwargs...)

    @eval begin function NLsolve.$(fn)(ms::SCPAProblem, ωs_init, ψs_init::VectorField; kwargs...)
            ms.m==size(ψs_init,2) || throw("number of modes in ψs_init ($(size(ψs_init,2))) must be the same as given in SCPA ($(ms.m))")
            return $(fn)(ms.SALT,ωs_init,ψs_init; kwargs...)
        end
    end
end


"""
    nlsolve(scpa; kwargs...) -> results

    nlsolve(scpa, init_ωs, init_ψs; kwargs...) -> results

solve nonlinear SCPA equation where `scpa` is a `HelmholtzSCPA` or `MaxwellSCPA`.
Results also stored `scpa`.

`init_ψs` must be a `ScalarField` (for Helmholtz) or `ElectricField` (for Maxwell)
"""
nlsolve

"""
    fixedpoint(scpa; kwargs...) -> results

    fixedpoint(scpa, init_ωs, init_ψs; kwargs...) -> results

Convenience wrapper for `NLsolve`'s `fixedpoint` wrapper.
Solves nonlinear SCPA equation where `scpa` is a `HelmholtzSCPA` or `MaxwellSCPA`.
Results also stored `scpa`.

`init_ψs` must be a `ScalarField` (for Helmholtz) or `ElectricField` (for Maxwell)
"""
fixedpoint

################################################################################
# Pretty Printing

import ..Common.PRINTED_COLOR_LIGHT

function Base.show(io::IO,ms::SCPAProblem)
    print(io,ms.m, " mode")
    ms.m>1 ? print(io,"s") : nothing
    if typeof(ms)<:HelmholtzSCPA
        printstyled(io," HelmholtzSCPA",color = PRINTED_COLOR_LIGHT)
    elseif typeof(ms)<:MAxwellSCPA
        printstyled(io," MaxwellSCPA",color = PRINTED_COLOR_LIGHT)
    end
    print(IOContext(io,:SCPA=>true), ms.SALT)
end

################################################################################
# Plotting
@recipe f(ms::SCPAProblem; by=abs2) = ms, by
@recipe f(by::Function, ms::SCPAProblem) = ms, by
@recipe f(ms::SCPAProblem, by::Function) = ms.SALT, by

@recipe f(sim::Simulation, ms::SCPAProblem; by=abs2) = sim, ms, by
@recipe f(ms::SCPAProblem, sim::Simulation; by=abs2) = sim, ms, by

@recipe f(ms::SCPAProblem, sim::Simulation, by::Function) = sim, ms, by
@recipe f(sim::Simulation, by::Function, ms::SCPAProblem) = sim, ms, by
@recipe f(ms::SCPAProblem, by::Function, sim::Simulation) = sim, ms, by
@recipe f(by::Function, sim::Simulation, ms::SCPAProblem) = sim, ms, by
@recipe f(by::Function, ms::SCPAProblem, sim::Simulation) = sim, ms, by

@recipe f(sim::Simulation, ms::SCPAProblem, by::Function) = sim, ms.SALT, by

end # module

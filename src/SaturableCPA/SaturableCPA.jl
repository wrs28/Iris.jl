# TODO: add simulation-modifying convenience wrappers to make bc's incoming. Without that, this module provides essentially no utility above SALT
"""
    module SaturableCPA

Convenience wrappers for solution of the nonlinear SatruableCPA equation.
Exports `MaxwellSCPA` constructor and extends methods of `NLsolve`, which must
be separately imported. Also convenience nonlinear solver wrapper `SCPA`.
"""
module SaturableCPA

export SCPA
export MaxwellSCPA

using ..Common
using ..Lasing
using NLsolve
using RecipesBase

struct MaxwellSCPA{TMS} SALT::TMS end

fnames = (:nlsolve,:fixedpoint)
for fn ∈ fnames
    @eval NLsolve.$(fn)(ms::MaxwellSCPA; kwargs...) = $(fn)(ms, ms.ωs, ms.ψs; kwargs...)

    @eval begin function NLsolve.$(fn)(ms::MaxwellSCPA, ωs_init, ψs_init::ElectricField; kwargs...)
            ms.m==size(ψs_init,2) || throw("number of modes in ψs_init ($(size(ψs_init,2))) must be the same as given in MaxwellSCPA ($(ms.m))")
            return $(fn)(ms.SALT,ωs_init,ψs_init; kwargs...)
        end
    end
end

function Base.getproperty(mcpa::MaxwellSCPA,sym::Symbol)
    if sym == :SALT
        getfield(mcpa,:SALT)
    else
        getproperty(getfield(mcpa,:SALT),sym)
    end
end

Base.propertynames(::MaxwellSCPA,private=false) = propertynames(MaxwellSALT,private)

"""
    SCPA(simulation, ωs_init, ψs_init; kwargs...) -> ωs, ψs

Convenience wrapper, solving SCPA problem for `simulation`, initialized at CPA
frequencies `ωs_init` and fields `ψs_init`.
"""
function SCPA(
            sim::Simulation,
            ωs_init::Vector{Float64},
            ψs_init::ElectricField;
            kwargs...)

    ms = MaxwellSCPA(sim,length(ωs_init))
    results = nlsolve(ms, ωs_init, ψs_init; kwargs...)
    ωs,ψs = x_to_ωψ(sim,results.zero,length(ωs_init))
    (results.f_converged || results.x_converged) || printstyled("beware: no convergence",color=PRINTED_COLOR_BAD)
    return ωs, ψs
end

################################################################################
# Pretty Printing

import ..Common.PRINTED_COLOR_LIGHT
import ..Common.PRINTED_COLOR_WARN
import ..Common.PRINTED_COLOR_VARIABLE
import ..Common.PRINTED_COLOR_GOOD
import ..Common.PRINTED_COLOR_BAD
import ..Common.PRINTED_COLOR_INSTRUCTION

function Base.show(io::IO,ms::MaxwellSCPA)
    print(io,ms.m, " mode")
    ms.m>1 ? print(io,"s") : nothing
    printstyled(io," MaxwellSCPA",color = PRINTED_COLOR_LIGHT)
    if !ms.solved[]
        printstyled(io," (",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"unsolved",color=PRINTED_COLOR_WARN)
        printstyled(io,", pass to ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"NLsolve.nlsolve",color=PRINTED_COLOR_VARIABLE)
        printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
    elseif ms.converged[]
        printstyled(io," (",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"solved",color=PRINTED_COLOR_GOOD)
        printstyled(io,", solution in fields ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"ωs",color=:cyan)
        printstyled(io,", ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"ψs",color=:cyan)
        printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
    else
        printstyled(io," (",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"solve attempted",color=PRINTED_COLOR_BAD)
        printstyled(io,", ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"no convergence",color=PRINTED_COLOR_BAD)
        printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
    end
end

# Plotting
@recipe f(ms::MaxwellSCPA;by=abs2) = ms,by
@recipe f(by::Function,ms::MaxwellSCPA) = ms,by
@recipe f(ms::MaxwellSCPA,by::Function) = ms.SALT,by
@recipe f(sim::Simulation,ms::MaxwellSCPA;by=abs2) = sim,ms.SALT,by
@recipe f(ms::MaxwellSCPA,sim::Simulation;by=abs2) = sim,ms.SALT,by

end # module

# TODO: add simulation-modifying convenience wrappers to make bc's incoming. Without that, this module provides essentially no utility above SALT
"""
    module SaturableCPA

Convenience wrappers for solution of the nonlinear SatruableCPA equation.
Exports `maxwell_scpa` constructor and extends methods of `NLsolve`, which must
be separately imported. Also convenience wrapper `SCPA`.
"""
module SaturableCPA

export SCPA
export maxwell_scpa

using ..Common
using ..Lasing
using NLsolve
using RecipesBase

import ..Common.PRINTED_COLOR_LIGHT
import ..Common.PRINTED_COLOR_WARN
import ..Common.PRINTED_COLOR_VARIABLE
import ..Common.PRINTED_COLOR_GOOD
import ..Common.PRINTED_COLOR_BAD
import ..Common.PRINTED_COLOR_INSTRUCTION

struct Maxwell_SCPA{TMS} SALT::TMS end

maxwell_scpa(args...; kwargs...) = Maxwell_SCPA(maxwell_salt(args...; kwargs...))


fnames = (:nlsolve,:fixedpoint)
for fn ∈ fnames

    @eval NLsolve.$(fn)(ms::Maxwell_SCPA; kwargs...) = $(fn)(ms, ms.ωs, ms.ψs; kwargs...)

    @eval begin function NLsolve.$(fn)(ms::Maxwell_SCPA, ωs_init::Union{Real,Vector}, ψs_init::ElectricField; kwargs...)
            N = ms.m
            n,m = size(ψs_init)
            N==m || throw("number of modes in ψs_init ($m) must be the same as given in Maxwell_SCPA ($N)")
            return $(fn)(ms.SALT,ωs_init,ψs_init; kwargs...)
        end
    end

end

function Base.getproperty(mcpa::Maxwell_SCPA,sym::Symbol)
    if sym == :SALT
        getfield(mcpa,:SALT)
    else
        getproperty(getfield(mcpa,:SALT),sym)
    end
end

Base.propertynames(::Maxwell_SCPA,private=false) = propertynames(Maxwell_SALT,private)

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

    ms = maxwell_scpa(sim,length(ωs_init))
    results = nlsolve(ms, ωs_init, ψs_init; kwargs...)
    ωs,ψs = x_to_ωψ(sim,results.zero,length(ωs_init))
    (results.f_converged || results.x_converged) || printstyled("beware: no convergence",color=PRINTED_COLOR_BAD)
    return ωs, ψs
end

# Pretty Printing
function Base.show(io::IO,ms::Maxwell_SCPA)
    print(io,ms.m, " mode")
    ms.m>1 ? print(io,"s") : nothing
    printstyled(io," Maxwell_SCPA",color = PRINTED_COLOR_LIGHT)
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
@recipe f(ms::Maxwell_SCPA;by=abs2) = ms,by
@recipe f(by::Function,ms::Maxwell_SCPA) = ms,by
@recipe f(ms::Maxwell_SCPA,by::Function) = ms.SALT,by
@recipe f(sim::Simulation,ms::Maxwell_SCPA;by=abs2) = sim,ms.SALT,by
@recipe f(ms::Maxwell_SCPA,sim::Simulation;by=abs2) = sim,ms.SALT,by

end # module

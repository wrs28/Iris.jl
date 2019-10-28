"""
    module SaturableCPA
"""
module SaturableCPA

export SCPA
export maxwell_scpa

using ..Common
using ..Lasing
using NLsolve
using RecipesBase

const PRINTED_COLOR_LIGHT = Common.PRINTED_COLOR_LIGHT
const PRINTED_COLOR_WARN = Common.PRINTED_COLOR_WARN
const PRINTED_COLOR_VARIABLE = Common.PRINTED_COLOR_VARIABLE
const PRINTED_COLOR_GOOD = Common.PRINTED_COLOR_GOOD
const PRINTED_COLOR_BAD = Common.PRINTED_COLOR_BAD

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

struct Maxwell_SCPA{TMS}
    SALT::TMS
end
maxwell_scpa(sim::Simulation,m::Int) = Maxwell_SCPA(maxwell_salt(conj(sim),m))


fnames = (:nlsolve,:fixedpoint)
for fn ∈ fnames
    @eval begin function NLsolve.$(fn)(ms::Maxwell_SCPA, ωs_init::Vector, ψs_init::ElectricField; kwargs...)
            return nlsolve(ms.SALT,ωs_init,conj(ψs_init);kwargs...)
        end
    end
end

function Base.getproperty(ms::Maxwell_SCPA,sym::Symbol)
    if sym == :ωs
        getfield(ms.SALT,:ωs)
    elseif sym == :ψs
        getfield(ms.SALT,:ψs)
    else
        return getfield(ms,sym)
    end
end

function Base.show(io::IO,ms::Maxwell_SCPA)
    printstyled(io,"Maxwell_SCPA ",color = PRINTED_COLOR_LIGHT)
    print(io,"with ", ms.SALT.m, " mode")
    ms.SALT.m>1 ? print(io,"s") : nothing
    if !ms.SALT.solved[]
        print(io," (")
        printstyled(io,"unsolved",color=PRINTED_COLOR_WARN)
        print(io,", for use in ")
        printstyled(io,"nlsolve",color=PRINTED_COLOR_VARIABLE)
        print(io,")")
    elseif ms.SALT.converged[]
        print(io," (")
        printstyled(io,"solved",color=PRINTED_COLOR_GOOD)
        print(io,", solution in fields ")
        printstyled(io,"ωs",color=:cyan)
        print(io,", ")
        printstyled(io,"ψs",color=:cyan)
        print(io,")")
    else
        print(io," (")
        printstyled(io,"solve attempted",color=PRINTED_COLOR_BAD)
        print(io,", ")
        printstyled(io,"no convergence",color=PRINTED_COLOR_BAD)
        print(io,")")
    end
end

# Plotting
@recipe f(ms::Maxwell_SCPA;by=abs2) = ms.ψs,by
@recipe f(sim::Simulation,ms::Maxwell_SCPA;by=abs2) = sim,ms.ψs,by
@recipe f(ms::Maxwell_SCPA,sim::Simulation;by=abs2) = sim,ms.ψs,by

end # module

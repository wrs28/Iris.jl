"""
    module SaturableCPA
"""
module SaturableCPA

export SCPA
export maxwell_scpa

using ..Common
using ..Lasing
using NLsolve

function SCPA(sim::Simulation, ωs_init::Array{Float64}, ψs_init::Array; kwargs...)
    csim = conj(sim)
    ωs,ψs = SALT(sim,ωs_init,conj(ψs_init);kwargs...)
    return ωs,conj(ψs)
end

struct Maxwell_SCPA{TMS}
    MS::TMS
end
maxwell_scpa(sim::Simulation,m::Int) = Maxwell_SCPA(maxwell_salt(conj(sim),m))


fnames = (:nlsolve,:fixedpoint)
for fn ∈ fnames
    @eval begin function NLsolve.$(fn)(ms::Maxwell_SCPA, ωs_init::Array, ψs_init::Array; kwargs...)
            return nlsolve(ms.MS,ωs_init,conj(ψs_init);kwargs...)
        end
    end
end

function Base.getproperty(ms::Maxwell_SCPA,sym::Symbol)
    if sym == :ωs
        getfield(ms.MS,:ωs)
    elseif sym == :ψs
        getfield(ms.MS,:ψs)
    else
        return getfield(ms,sym)
    end
end

function Base.show(io::IO,mscpa::Maxwell_SCPA)
    ms = mscpa.MS
    print(io,"Maxwell_SCPA with $(ms.m) mode")
    ms.m>1 ? print(io,"s") : nothing
    if !ms.solved[]
        print(io," (")
        printstyled(io,"unsolved",color=:light_yellow)
        print(io,", for use in ")
        printstyled(io,"nlsolve",color=:cyan)
        print(io,")")
    elseif ms.converged[]
        print(io," (")
        printstyled(io,"solved",color=:light_green)
        print(io,", solution in fields ")
        printstyled(io,"ωs",color=:cyan)
        print(io,", ")
        printstyled(io,"ψs",color=:cyan)
        print(io,")")
    else
        print(io," (")
        printstyled(io,"solve attempted",color=:light_red)
        print(io,", ")
        printstyled(io,"no convergence",color=:light_red)
        print(io,")")
    end
end

end # module

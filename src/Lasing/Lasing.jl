module Lasing

export SALT
export maxwell_salt

const INDEX_OFFSET = 7

using ..Common
using NLsolve


"""
    SALT(simulation, ωs_init, ψs_init; kwargs...) -> ωs,ψs
"""
SALT

"""
    Maxwell_SALT(simulation, m::Int) -> ms
"""
Maxwell_SALT

function SALT(sim::Simulation, ωs_init::Array{Float64}, ψs_init::Array; kwargs...)
    results = nlsolve(sim, ωs_init, ψs_init; kwargs...)
    ωs,ψs = x_to_ωψ(results.zero,length(ωs_init))
    (results.f_converged || results.x_converged) || printstyled("beware: no convergence",color=:red)
    return ωs, ψs
end

function SALT_bootstrap(sim::Simulation,D0)

struct Maxwell_SALT{TM}
    M::TM
    ωs::Array{Float64,1}
    ψs::Array{ComplexF64,2}
    x::Array{Float64,1}
    n::Int
    m::Int
    res::Array{ComplexF64,2}
    index::Int
    solved::Ref{Bool}
    converged::Ref{Bool}
end
function Maxwell_SALT(sim::Simulation,m::Int)
    M = maxwell(sim)
    n = 3length(sim.x)
    ωs = Array{Float64,1}(undef,m)
    ψs = Array{ComplexF64,2}(undef,n,m)
    x = ωψ_to_x(ωs,ψs)
    res = Array{ComplexF64,2}(undef,n,m)
    index = n÷2+INDEX_OFFSET
    return Maxwell_SALT(M,ωs,ψs,x,n,m,res,index,Ref(false),Ref(false))
end
maxwell_salt(sim::Simulation,m::Int) = Maxwell_SALT(sim,m)

fnames = (:nlsolve,:fixedpoint)
for fn ∈ fnames
    @eval NLsolve.$(fn)(sim::Simulation, ωs_init::Array, ψs_init::Array; kwargs...) = $(fn)(Maxwell_SALT(sim,length(ωs_init)), ωs_init, ψs_init; kwargs...)
    @eval begin function NLsolve.$(fn)(ms::Maxwell_SALT, ωs_init::Array, ψs_init::Array; kwargs...)
            ωψ_to_x!(ms.x,ωs_init,ψs_init)
            results =  $(fn)(ms,ms.x;kwargs...)
            x_to_ωψ!(ms,results.zero)
            ms.solved[] = true
            ms.converged[] = results.f_converged || results.x_converged
            return results
        end
    end
end

@inline function (ms::Maxwell_SALT)(F,x)
    x_to_ωψ!(ms,x)
    for μ ∈ eachindex(ms.ωs)
        ms.res[:,μ] = ms.M(ms.ωs[μ],ms.ωs,ms.ψs)*ms.ψs[:,μ]/ms.ψs[ms.index,μ]
    end
    for i ∈ eachindex(ms.res)
        F[i], F[ms.n*ms.m + i] = reim(ms.res[i])
    end
    return nothing
end

function Base.show(io::IO,ms::Maxwell_SALT)
    print(io,"Maxwell_SALT with $(ms.m) mode")
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

################################################################################

function ωψ_to_x(ωs,ψs)
    n,m = size(ψs)
    x = Array{Float64,1}(undef,2m*n)
    ωψ_to_x!(x,ωs,ψs)
    return x
end
ωψ_to_x!(ms::Maxwell_SALT,ωs,ψs) = ωψ_to_x!(ms.x,ωs,ψs)
@inline function ωψ_to_x!(x,ωs,ψs)
    n,m = size(ψs)
    nm = n*m
    index = n÷2+INDEX_OFFSET
    for μ ∈ 1:m
        for i ∈ 1:n
            if i == index
                x[(μ-1)n+index] = abs(ψs[index,μ])
                x[nm+(μ-1)n+index] = ωs[μ]
            else
                x[(μ-1)n+i], x[nm+(μ-1)n+i] = reim(ψs[i,μ]/ψs[index,μ])
            end
        end
    end
    return nothing
end


function x_to_ωψ(x,m)
    n = length(x)÷2÷m
    ωs = Array{Float64,1}(undef,m)
    ψs = Array{ComplexF64,2}(undef,n,m)
    x_to_ωψ!(ωs,ψs,x)
    return ωs,ψs
end
x_to_ωψ!(ms::Maxwell_SALT,x) = x_to_ωψ!(ms.ωs,ms.ψs,x)
@inline function x_to_ωψ!(ωs,ψs,x)
    n,m = size(ψs)
    nm = n*m
    index = n÷2+INDEX_OFFSET
    for μ ∈ 1:m
        for i ∈ 1:n
            if i == index
                ψs[index,μ] = complex(x[(μ-1)n+index],0)
            else
                ψs[i,μ] = complex(x[(μ-1)n+i],x[nm+(μ-1)n+i])*x[(μ-1)n+index]
            end
        end
        ωs[μ] = x[nm + (μ-1)n+index]
    end
    return nothing
end


end # module

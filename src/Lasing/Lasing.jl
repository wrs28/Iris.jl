module Lasing

export SALT
export maxwell_salt

const INDEX_OFFSET = 7

using ..Common
using LinearAlgebra
using NLsolve
using RecipesBase

const PRINTED_COLOR_LIGHT = Common.PRINTED_COLOR_LIGHT
const PRINTED_COLOR_WARN = Common.PRINTED_COLOR_WARN
const PRINTED_COLOR_VARIABLE = Common.PRINTED_COLOR_VARIABLE
const PRINTED_COLOR_GOOD = Common.PRINTED_COLOR_GOOD
const PRINTED_COLOR_BAD = Common.PRINTED_COLOR_BAD

"""
    SALT(simulation, ωs_init, ψs_init; kwargs...) -> ωs,ψs
"""
SALT

"""
    Maxwell_SALT(simulation, m::Int) -> ms
"""
Maxwell_SALT

function SALT(
            sim::Simulation,
            ωs_init::Vector{Float64},
            ψs_init::ElectricField;
            kwargs...)

    ms = maxwell_salt(sim,length(ωs_init))
    results = nlsolve(ms, ωs_init, ψs_init; kwargs...)
    ωs,ψs = x_to_ωψ(sim,results.zero,length(ωs_init))
    (results.f_converged || results.x_converged) || printstyled("beware: no convergence",color=PRINTED_COLOR_BAD)
    return ωs, ψs
end

function SALT_bootstrap(sim::Simulation,D0)

end

struct Maxwell_SALT{NMODES,TM,DIM}
    M::TM
    ωs::Vector{Float64}
    ψs::ElectricField{DIM}
    x::Vector{Float64}
    n::Int
    m::Int
    res::Matrix{ComplexF64}
    index::Int
    solved::Ref{Bool}
    converged::Ref{Bool}
end
function Maxwell_SALT(sim::Simulation{DIM},m::Int) where DIM
    M = maxwell(sim)
    n = 3length(sim)
    ωs = Vector{Float64}(undef,m)
    ψs = ElectricField(sim.x,m)
    x = ωψ_to_x(ωs,ψs)
    res = Matrix{ComplexF64}(undef,n,m)
    index = n÷2+INDEX_OFFSET
    return Maxwell_SALT{m,typeof(M),DIM}(M,ωs,ψs,x,n,m,res,index,Ref(false),Ref(false))
end
maxwell_salt(sim::Simulation,m::Int) = Maxwell_SALT(sim,m)

fnames = (:nlsolve,:fixedpoint)
for fn ∈ fnames
    # @eval NLsolve.$(fn)(sim::Simulation, ωs_init::Vector, ψs_init::ElectricField; kwargs...) = $(fn)(Maxwell_SALT(sim,length(ωs_init)), ωs_init, ψs_init; kwargs...)
    @eval begin function NLsolve.$(fn)(ms::Maxwell_SALT{N}, ωs_init::Vector, ψs_init::ElectricField; kwargs...) where N
            n,m = size(ψs_init)
            N==m || throw("number of modes in initial ψs ($n) must be the same as given in Maxwell_SALT ($N)")
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
        res = view(ms.res,:,μ)
        mul!(res,ms.M(ms.ωs[μ],ms.ωs,ms.ψs),ms.ψs[:,μ])
        rdiv!(res,ms.ψs[ms.index,μ])
    end
    for i ∈ eachindex(ms.res) F[i], F[ms.n*ms.m + i] = reim(ms.res[i]) end
    return nothing
end
# @inline function (ml::Maxwell_NLS)(F,J,x)
#     x_to_ψ!(ml,x)
#     mul!(ml.res,ml.M(ml.ω,[ml.ω],ml.ψ),ml.ψ.val)
#     if !isnothing(F)
#         for i ∈ eachindex(ml.res) F[i], F[length(x)÷2 + i] = reim(ml.res[i] - ml.j[i]) end
#     end
#     if !isnothing(J)
#         anchor = copy(ml.res)
#         temp1 = similar(anchor)
#         temp2 = similar(anchor)
#         for j ∈ 1:length(temp1)
#             ml.ψ.val[j] += 1e-5
#             mul!(temp1,ml.M(ml.ω,[ml.ω],ml.ψ),ml.ψ.val)
#             ml.ψ.val[j] -= 1e-5
#             temp2[:] = (temp1-anchor)/1e-5
#             for i ∈ eachindex(anchor) J[i,j], J[length(x)÷2 + i,j] = reim(temp2[i]) end
#
#             ml.ψ.val[j] += 1e-5im
#             mul!(temp1,ml.M(ml.ω,[ml.ω],ml.ψ),ml.ψ.val)
#             ml.ψ.val[j] -= 1e-5im
#             temp2[:] = (temp1-anchor)/1e-5
#             for i ∈ eachindex(anchor) J[i,length(x)÷2 + j], J[length(x)÷2 + i,length(x)÷2 + j] = reim(temp2[i]) end
#         end
#     end
#     return nothing
# end

function Base.show(io::IO,ms::Maxwell_SALT)
    printstyled(io,"Maxwell_SALT ",color = PRINTED_COLOR_LIGHT)
    print(io,"with ", ms.m, " mode")
    ms.m>1 ? print(io,"s") : nothing
    if !ms.solved[]
        print(io," (")
        printstyled(io,"unsolved",color=PRINTED_COLOR_WARN)
        print(io,", for use in ")
        printstyled(io,"nlsolve",color=PRINTED_COLOR_VARIABLE)
        print(io,")")
    elseif ms.converged[]
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

################################################################################

function ωψ_to_x(ωs::Vector,ψs::ElectricField)
    n,m = size(ψs)
    x = Vector{Float64}(undef,2m*n)
    ωψ_to_x!(x,ωs,ψs)
    return x
end
ωψ_to_x!(ms::Maxwell_SALT,ωs,ψs) = ωψ_to_x!(ms.x,ωs,ψs)
ωψ_to_x!(x,ωs::Vector,ψs::ElectricField) = ωψ_to_x!(x,ωs,ψs.values)
@inline function ωψ_to_x!(x,ωs::Vector,ψs::Matrix)
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


function x_to_ωψ(sim::Simulation{N},x,m) where N
    ωs = Vector{Float64}(undef,m)
    ψs = ElectricField(sim.x,m)
    x_to_ωψ!(ωs,ψs,x)
    return ωs,ψs
end
x_to_ωψ!(ms::Maxwell_SALT,x) = x_to_ωψ!(ms.ωs,ms.ψs,x)
@inline function x_to_ωψ!(ωs::Vector,ψs::ElectricField,x)
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

# Plotting
@recipe f(ms::Maxwell_SALT;by=abs2) = ms.ψs,by
@recipe f(sim::Simulation,ms::Maxwell_SALT;by=abs2) = sim,ms.ψs,by
@recipe f(ms::Maxwell_SALT,sim::Simulation;by=abs2) = sim,ms.ψs,by

end # module

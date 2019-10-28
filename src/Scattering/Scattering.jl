"""
    module Scattering
"""
module Scattering

files = (
    "1D/Scattering1D.jl",
    # "2D/Scattering2D.jl",
    # "3D/Scattering3D.jl",
    )

export scattering
export scattering_nl
export maxwell_nls

using ..Common
using ..Spectral
using LinearAlgebra
using NLsolve
using SparseArrays

const EQUIVALENT_SOURCE_RELATIVE_CUTOFF = 1e-8

const PRINTED_COLOR_LIGHT = Common.PRINTED_COLOR_LIGHT
const PRINTED_COLOR_NUMBER = Common.PRINTED_COLOR_NUMBER
const PRINTED_COLOR_BAD = Common.PRINTED_COLOR_BAD
const PRINTED_COLOR_WARN = Common.PRINTED_COLOR_WARN
const PRINTED_COLOR_VARIABLE = Common.PRINTED_COLOR_VARIABLE

struct ScatteringSolution{N}
    total::ElectricField{N}
    incident::ElectricField{N}
    scattered::ElectricField{N}
    ω::Float64
    ScatteringSolution(t::T,i::T,s::T,ω) where T<:ElectricField{N} where N = new{N}(t,i,s,ω)
end
ScatteringSolution(pos::Array{T,1},tot::Array,inc::Array,sct::Array,ω::Real) where T<:Point =
    ScatteringSolution(ElectricField(pos,tot),ElectricField(pos,inc),ElectricField(pos,sct),ω)

struct EquivalentSource{N,TS,TSIM}
    mask::TS
    field::ElectricField{N}
    simulation::TSIM
    EquivalentSource(mask::AbstractShape{N},field,sim::TSIM) where {N,TSIM<:Simulation} = new{N,typeof(mask),TSIM}(mask,field,sim)
end

struct Source{N}
    p::Point{N}
    j::SparseVector{ComplexF64,Int}
    Source(p::Point{N},j) where N = new{N}(p,j)
end

struct Maxwell_NLS{N,TM,TS}
    ω::Float64
    M::TM
    equivalent_source::TS
    j::SparseMatrixCSC{ComplexF64,Int}
    ψ::ElectricField{N}
    residual::Matrix{ComplexF64}
    x::Vector{Float64}
    solved::Ref{Bool}
    converged::Ref{Bool}
end
@inline function (ml::Maxwell_NLS)(F,J,x)
    x_to_ψ!(ml,x)
    A = ml.M(ml.ω,[ml.ω],ml.ψ)
    mul!(ml.res,A,ml.ψ.val)
    if !isnothing(F)
        for i ∈ eachindex(ml.res) F[i], F[length(x)÷2 + i] = reim(ml.res[i] - ml.j[i]) end
    end
    if !isnothing(J)
        anchor = copy(ml.res)
        temp1 = similar(anchor)
        temp2 = similar(anchor)
        for j ∈ 1:length(temp1)
            ml.ψ.val[j] += 1e-7
            mul!(temp1,ml.M(ml.ω,[ml.ω],ml.ψ),ml.ψ.val)
            ml.ψ.val[j] -= 1e-7
            temp2[:] = (temp1-anchor)/1e-7
            for i ∈ eachindex(anchor) J[i,j], J[length(x)÷2 + i,j] = reim(temp2[i]) end

            ml.ψ.val[j] += 1e-7im
            mul!(temp1,ml.M(ml.ω,[ml.ω],ml.ψ),ml.ψ.val)
            ml.ψ.val[j] -= 1e-7im
            temp2[:] = (temp1-anchor)/1e-7
            for i ∈ eachindex(anchor) J[i,length(x)÷2 + j], J[length(x)÷2 + i,length(x)÷2 + j] = reim(temp2[i]) end
        end
    end
    return nothing
end
maxwell_nls(sim::Simulation,args...;kwargs...) = Maxwell_NLS(sim,args...; kwargs...)

function NLsolve.nlsolve(nls::Maxwell_NLS{N},ψ_init::ElectricField{N};kwargs...) where N
    x_init = ψ_to_x(ψ_init)
    return nlsolve(only_fj!(nls),x_init;kwargs...)
end

################################################################################

function Base.getproperty(sct::ScatteringSolution,sym::Symbol)
    if Base.sym_in(sym,(:incoming,:inc,:Incident,:Inc,:In,:IN,:in,:i,:I))
        return getfield(sct,:incident)
    elseif Base.sym_in(sym,(:Total,:tot,:Tot,:t,:T))
        return getfield(sct,:total)
    elseif Base.sym_in(sym,(:Scattered,:scat,:Scat,:sct,:SCT,:Sct,:S,:s,:outgoing,:out,:Out,:OUT,:o,:O))
        return getfield(sct,:scattered)
    elseif Base.sym_in(sym,(:frequency,:freq,:omega,:Omega,:Freq,:Frequency,:f,:F))
        return getfield(sct,:ω)
    else
        return getfield(sct,sym)
    end
end

function Base.getproperty(es::EquivalentSource,sym::Symbol)
    if sym==:sim
        return getfield(es,:simulation)
    else
        return getfield(es,sym)
    end
end

function Base.getproperty(es::Maxwell_NLS,sym::Symbol)
    if sym==:res
        return getfield(es,:residual)
    else
        return getfield(es,sym)
    end
end

################################################################################

foreach(include,files)

################################################################################

function ψ_to_x(ψ::ElectricField)
    x = Vector{Float64}(undef,2prod(size(ψ)))
    ψ_to_x!(x,ψ)
    return x
end
ψ_to_x!(ml::Maxwell_NLS,ψ) = ψ_to_x!(ml.x,ψ)
ψ_to_x!(x,ψ::ElectricField) = ψ_to_x!(x,ψ.values)
@inline function ψ_to_x!(x,ψ::Matrix)
    n,m = size(ψ)
    nm = n*m
    for μ ∈ 1:m for i ∈ 1:n
        x[(μ-1)n+i], x[nm+(μ-1)n+i] = reim(ψ[i,μ])
    end end
    return nothing
end

function x_to_ψ(sim::Simulation,x,m)
    ψ = ElectricField(sim.x,m)
    x_to_ψ!(ψ,x)
    return ψ
end
x_to_ψ!(ml::Maxwell_NLS,x) = x_to_ψ!(ml.ψ,x)
@inline function x_to_ψ!(ψ::ElectricField,x)
    n,m = size(ψ)
    nm = n*m
    for μ ∈ 1:m for i ∈ 1:n
            ψ[i,μ] = complex(x[(μ-1)n+i],x[nm+(μ-1)n+i])
    end end
    return nothing
end

################################################################################
# Pretty Printing

function Base.show(io::IO,sc::EquivalentSource{N}) where N
    print(io,"$(N)D ")
    printstyled(io,"EquivalentSource",color=PRINTED_COLOR_LIGHT)
    print(io," (call with argument ")
    printstyled(io,"ω",color=PRINTED_COLOR_NUMBER)
    print(io,")")
end

function Base.show(io::IO,s::Source{N}) where N
    print(io,"$(N)D ")
    printstyled(io,"Source",color=PRINTED_COLOR_LIGHT)
    print(io," at frequency ")
    printstyled(io,"ω=",sc.ω,color=PRINTED_COLOR_NUMBER)
end

function Base.show(io::IO,sc::ScatteringSolution{N}) where N
    print(io,"$(N)D ")
    printstyled(io,"ScatteringSolution",color=PRINTED_COLOR_LIGHT)
    print(io," at frequency ")
    printstyled(io,"ω=",sc.ω,color=PRINTED_COLOR_NUMBER)
end

function Base.show(io::IO,ms::Maxwell_NLS)
    printstyled(io,"Maxwell_NLS ",color = PRINTED_COLOR_LIGHT)
    print(io,"nonlinear scattering problem")
    if !ms.solved[]
        print(io," (")
        printstyled(io,"unsolved",color=PRINTED_COLOR_WARN)
        print(io,", for use in ")
        printstyled(io,"nlsolve",color=PRINTED_COLOR_VARIABLE)
        print(io,")")
    elseif ms.converged[]
        print(io," (")
        printstyled(io,"solved",color=PRINTED_COLOR_GOOD)
        print(io,", solution in field ")
        printstyled(io,"ψ",color=:cyan)
        print(io,")")
    else
        print(io," (")
        printstyled(io,"solve attempted",color=PRINTED_COLOR_BAD)
        print(io,", ")
        printstyled(io,"no convergence",color=PRINTED_COLOR_BAD)
        print(io,")")
    end
end

end # module

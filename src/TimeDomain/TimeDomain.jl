"""
    module TimeDomain
"""
module TimeDomain

files = (
         "1D/TimeDomain1D.jl",
         # "2D/TimeDomain2D.jl",
         # "3D/TimeDomain3D.jl",
        )

export HelmholtzFDTD
export HelmholtzBlochFDTD
export propagate!

using ..Common
using Formatting
using RecipesBase
import ..Common.DEFAULT_CFL_NUMBER
import LinearAlgebra: mul!

abstract type AbstractFields{N,M} end
abstract type AbstractFDTD{N,M} end

include("DivGradCurls.jl")

propagate!() = nothing

################################################################################
# Define Field Structures

struct HelmholtzWaveFields{N,M} <: AbstractFields{N,M} # N is dimension, M is memory depth
    φ::NTuple{M,VectorField{N,1,Float64}}
    ∇Φ::NTuple{N,NTuple{M,VectorField{N,1,Float64}}}
end

# struct HelmholtzMatterFields{N} <: AbstractFields{N,2}
#     ρ::VectorField{N,1,Float64}
#     ι::VectorField{N,1,Float64}
#     D::VectorField{N,1,Float64}
# end
#
# function Base.getproperty(hmf::HelmholtzMatterFields, sym::Symbol)
#     if sym==:P
#         return complex.(getfield(hmf,:ρ),getfield(hmf,:ι))
#     else
#         return getfield(hmf,sym)
#     end
# end
#
# function Base.propertynames(::HelmholtzMatterFields, private=false)
#     if private
#         return (:P, fieldnames(HelmholtzMatterFields)...)
#     else
#         return (:P,:D)
#     end
# end


# struct HelmholtzBlochFields{N} <: AbstractFields{N}
#     wave::HelmholtzWaveFields{N}
#     matter::HelmholtzMatterFields{N}
# end
#
# function Base.getproperty(hbf::HelmholtzBlochFields, sym::Symbol)
#     if Base.sym_in(sym, propertynames(getfield(hbf,:matter)), true)
#         return getproperty(getfield(hbf,:matter),sym)
#     elseif Base.sym_in(sym, propertynames(getfield(hbf,:wave)))
#         return getproperty(getfield(hbf,:wave),sym)
#     else
#         return getfield(hbf,sym)
#     end
# end
#
# function Base.propertynames(hbf::HelmholtzBlochFields, private=false)
#     if private
#         return (fieldnames(HelmholtzBlochField)..., propertynames(hbf.wave)..., propertynames(hbf.matter))
#     else
#         return (propertynames(hbf.wave)..., propertynames(hbf.matter))
#     end
# end

################################################################################
# HelmholtzFDTD

struct HelmholtzFDTD{N,M,TS} <: AbstractFDTD{N,M} # N is dimension, M is memory depth
    fields::HelmholtzWaveFields{N,M}
    grad::Gradient{N}
    div::Divergence{N}
    n::Vector{Int}
    t::Vector{Float64}
    dt::Float64
    source::TS
    α::Vector{Vector{Float64}}
    β::Vector{Vector{Float64}}
end

_getnumsteps(::AbstractFields{N,M}) where {N,M} = M

"""
    HelmholtzFDTD(sim; dt=sim.dx*$DEFAULT_CFL_NUMBER)
"""
function HelmholtzFDTD(sim::Simulation{N,Common.Symmetric}, dt::Real) where N
    fields = HelmholtzWaveFields(sim)
    grad = Gradient(sim)
    div = Divergence(sim)
    n = Int[0]
    t = Float64[dt/2]
    sc = 0
    α, β = _αβ(sim, dt)
    return HelmholtzFDTD{N,_getnumsteps(fields),typeof(sc)}(fields, grad, div, n, t, dt, sc, α, β)
end

################################################################################
# HelmholtzBlochFDTD
#
# struct HelmholtzBlochFDTD{N,TS} <: AbstractFDTD{N}
#     fields::HelmholtzBlochField{N}
#     differentials::HelmholtzBlochField{N}
#     grad::Gradient{N}
#     div::Divergence{N}
#     n::Vector{Int}
#     t::Vector{Float64}
#     dt::Float64
#     source::TS
#     ρE::Vector{Float64}
#     σEp1⁻¹dt::Vector{Float64}
#     ρG::Vector{Float64}
#     σGp1⁻¹dt::Vector{Float64}
#     ε⁻¹σEp1⁻¹dt::Vector{Float64}
#     Ω::Matrix{Float64}
# end
#
# """
#     HelmholtzBlochFDTD(sim; dt=sim.dx*$DEFAULT_CFL_NUMBER)
# """
# function HelmholtzBlochFDTD(sim::Simulation{N,Common.Symmetric}; dt::Real=sim.dx*DEFAULT_CFL_NUMBER) where N
#     fields = HelmholtzField(sim)
#     differentials = HelmholtzField(sim)
#     grad = Gradient(sim)
#     div = Divergence(sim)
#     n = Int[0]
#     t = Float64[dt/2]
#     sc = 0
#
#     σ = real(sim.σ[1])
#     σEp1⁻¹dt = dt ./(1 .+ σ*dt/2)
#     ρE = (1 .- σ*dt/2).*σEp1⁻¹
#
#     σ = real(sim.σ_half[1])
#     σGp1⁻¹dt = dt ./(1 .+ σ*dt/2)
#     ρG = (1 .- σ*dt/2).*σGp1⁻¹
#
#     return HelmholtzBlochFDTD{Ns,typeof(sc)}(fields, differentials, grad, div, n, t, dt, sc, ρE,
#                 σEp1⁻¹dt, ρG, σGp1⁻¹dt, σEp1⁻¹dt./sim.ε)
# end
#
#
# # @inline
# function propagate!(fdtd::HelmholtzBlochFDTD, n::Integer=1)
#     dt = fdtd.dt
#
#     E = fdtd.fields.E
#     dE = fdtd.differentials.E
#     ρE = fdtd.ρE
#     grad = fdtd.grad.components
#
#     G = fdtd.fields.G
#     dG = fdtd.differentials.G
#     ρG = fdtd.ρG
#     div = fdtd.div
#
#     ε⁻¹σEp1⁻¹dt = fdtd.ε⁻¹σEp1⁻¹dt
#     σGp1⁻¹dt = fdtd.σGp1⁻¹dt
#
#     foreach(1:n) do i
#         # @fastmath
#         mul!(dE, div, G)
#         # @fastmath @inbounds @simd
#         for j ∈ eachindex(dE) E[j] = ρE[j]*E[j] + ε⁻¹σEp1⁻¹dt[j]*dE[j] end
#         # @fastmath @inbounds
#         fdtd.t[1] += dt/2
#
#         # @fastmath
#         mul!(dG, grad, E)
#         # @fastmath @inbounds @simd
#         for j ∈ eachindex(dG) G[j] = ρG[j]*G[j] + σGp1⁻¹dt[j]*dG[j] end
#         # @fastmath @inbounds
#         fdtd.t[1] += dt/2
#     end
#     fdtd.n[1] += n
#     return nothing
# end

################################################################################

Base.eltype(af::AbstractFields) = Float64
Base.length(af::AbstractFields) = length(af.φ[1])
Base.ndims(af::AbstractFields{N}) where N = N

Base.eltype(fdtd::AbstractFDTD) = Float64
Base.length(fdtd::AbstractFDTD) = length(fdtd.fields)
Base.ndims(fdtd::AbstractFDTD{N}) where N = N

################################################################################
# Pretty Printing

import ..Common.PRINTED_COLOR_DARK
import ..Common.PRINTED_COLOR_LIGHT
import ..Common.PRINTED_COLOR_NUMBER

function Base.show(io::IO,af::AbstractFields{N}) where N
    if typeof(af)<:HelmholtzWaveFields
        printstyled(io,"HelmholtzFields",color=PRINTED_COLOR_DARK)
    elseif typeof(af)<:HelmholtzBlochFields
        printstyled(io,"HelmholtzBlochFields",color=PRINTED_COLOR_DARK)
    end
    print(io," (")
    printstyled(io, length(af), color=PRINTED_COLOR_NUMBER)
    print(io," points)")
end

function Base.show(io::IO,fdtd::AbstractFDTD{N}) where N
    print(io,"$(N)D ")
    if typeof(fdtd)<:HelmholtzFDTD
        printstyled(io,"HelmholtzFDTD",color=PRINTED_COLOR_LIGHT)
    end
    print(io," (")
    printstyled(io, length(fdtd), color=PRINTED_COLOR_NUMBER)
    print(io," points, step ")
    printstyled(io,"n=$(fdtd.n[1])", color=PRINTED_COLOR_NUMBER)
    print(io,", time ")
    printstyled(io,"t=",fmt("1.2e",fdtd.t[1]), color=PRINTED_COLOR_NUMBER)
    print(io,")")
end

################################################################################
# Plotting

@recipe function f(field::AbstractFields, by::Function)
    @series begin field.φ[1], by end
end

@recipe function f(fdtd::AbstractFDTD, by::Function)
    @series begin fdtd.fields, by end
end

@recipe f(v::AbstractFields; by=abs2) = v, by
@recipe f(by::Function,e::AbstractFields) = e,by

@recipe f(v::AbstractFDTD; by=abs2) = v, by
@recipe f(by::Function,e::AbstractFDTD) = e,by

################################################################################

foreach(include,files)

end # module

"""
    module TimeDomain
"""
module TimeDomain

files = (
         "1D/TimeDomain1D.jl",
         "2D/TimeDomain2D.jl",
         # "3D/TimeDomain3D.jl",
        )

export Constant
export HelmholtzFDTD
export HelmholtzBlochFDTD
export HelmholtzPointSource
export propagate!
export initialize!

using ..Common
using Formatting
using Plots
using RecipesBase

import ..Common.DEFAULT_CFL_NUMBER
import LinearAlgebra: mul!

abstract type AbstractFields{N,M} end
abstract type AbstractFDTD{N,M,P} end

include("DivGradCurls.jl")

propagate!() = nothing
_αβ() = nothing

################################################################################
# Define Field Structures

struct HelmholtzWaveFields{N,M} <: AbstractFields{N,M} # N is dimension, M is memory depth
    φ::NTuple{M,VectorField{N,1,Float64}}
    ∇Φ::NTuple{M,NTuple{N,VectorField{N,1,Float64}}}
end

function Base.getproperty(hwf::HelmholtzWaveFields{N,2}, sym::Symbol) where N
    if Base.sym_in(sym,(:φ,:ϕ,:phi,:psi,:ψ))
        return getfield(hwf,sym)[1]
    else
        return getfield(hwf,sym)
    end
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

struct AnimationOptions{TB,TO}
    animation::Animation
    interval::Int
    start::Int
    by::TB
    ylims::Vector{Float64}
    plotoptions::TO
end

function AnimationOptions(;
            interval::Integer=100,
            start::Integer=1,
            by::Function=real,
            ylims::Vector=Float64[0,0],
            plotoptions...)
    AnimationOptions(Animation(), interval, start, by, ylims, plotoptions)
end

struct HelmholtzFDTD{N,M,P,TS,TSIM} <: AbstractFDTD{N,M,P} # N is dimension, M is memory depth
    fields::HelmholtzWaveFields{N,M}
    grad::Gradient{N}
    div::Divergence{N}
    n::Vector{Int}
    t::Vector{Float64}
    dt::Float64
    source::TS
    α::Vector{Vector{Float64}}
    β::Vector{Vector{Float64}}
    pmldims::NTuple{P,Int}
    options::AnimationOptions
    simulation::TSIM
end

function HelmholtzFDTD(
            sim::Simulation{N,Common.Symmetric},
            sc::TS,
            dt::Real;
            kwargs...
            ) where {N,TS}

    fields = HelmholtzWaveFields(sim)
    grad = Gradient(sim)
    div = Divergence(sim)
    n = Int[0]
    t = Float64[dt/2]
    α, β = _αβ(sim, dt)
    pmldims = _whichdimspml(sim)
    return HelmholtzFDTD{N,_memorydepth(fields),_ndimpml(sim),TS,typeof(sim)}(fields, grad, div, n, t, dt, sc, α, β, pmldims, AnimationOptions(kwargs...), sim)
end

function Base.getproperty(fdtd::HelmholtzFDTD, sym::Symbol)
    if Base.sym_in(sym,(:φ,:ϕ,:phi,:ψ,:psi))
        return getproperty(getfield(fdtd,:fields),:φ)
    else
        return getfield(fdtd,sym)
    end
end

function Base.propertynames(::HelmholtzFDTD, private=false)
    if private
        return (:φ, fieldnames(HelmholtzFDTD)...)
    else
        return (:φ, :dt, :n, :options, :t)
    end
end


struct HelmholtzPointSource{FX,FW,FA,FP}
    xoft::FX    # position of time
    aoft::FA    # amplitude of time
    ωoft::FW    # frequency of time
    ϕoft::FP    # phase of time
    σ²::Float64
    N::Float64

    function HelmholtzPointSource(xoft, aoft, ωoft, ϕoft, σ²::Real, N::Real)
        typeof(xoft)<:Real ? xoft = Constant(xoft) : nothing
        typeof(aoft)<:Real ? aoft = Constant(aoft) : nothing
        typeof(ωoft)<:Real ? ωoft = Constant(ωoft) : nothing
        typeof(ϕoft)<:Real ? ϕoft = Constant(ϕoft) : nothing
        return new{typeof(xoft),typeof(aoft),typeof(ωoft),typeof(ϕoft)}(xoft,aoft,ωoft,ϕoft,σ²,N)
    end
end

function (ps::HelmholtzPointSource)(x::Point,t::Real)
    x0 = ps.xoft(t)
    typeof(x0)<:Point ? nothing : x0 = Point(x0...)
    ω = ps.ωoft(t)
    a = ps.aoft(t)
    ϕ = ps.ϕoft(t)
    r² = (x-x0).r^2
    return -a*sin(ϕ + ω*t) * exp(-r²/2ps.σ²)/ps.N/ω
end

function Base.show(io::IO, mime::MIME"text/plain", ps::HelmholtzPointSource)
    printstyled(io,"HelmholtzPointSource: "; color=PRINTED_COLOR_DARK)
    println(io)
    printstyled(io,"\tx(t)\t"; color=PRINTED_COLOR_VARIABLE)
    show(io,mime,ps.xoft)
    printstyled(io,"\n\tω(t)\t"; color=PRINTED_COLOR_VARIABLE)
    show(io,mime,ps.ωoft)
    printstyled(io,"\n\ta(t)\t"; color=PRINTED_COLOR_VARIABLE)
    show(io,mime,ps.aoft)
    printstyled(io,"\n\tϕ(t)\t"; color=PRINTED_COLOR_VARIABLE)
    show(io,mime,ps.ϕoft)
end

struct Constant<:Function c::Float64 end

(C::Constant)(t::Real) = C.c

Base.iszero(ps::HelmholtzPointSource) = iszero(ps.aoft)
Base.iszero(C::Constant) = iszero(C.c)

function Base.show(io::IO, ::MIME"text/plain", C::Constant)
    printstyled(io,C.c; color=PRINTED_COLOR_NUMBER)
    print(io," (constant)")
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
Base.length(af::AbstractFields) = length(af.φ)
Base.ndims(af::AbstractFields{N}) where N = N

Base.eltype(fdtd::AbstractFDTD) = Float64
Base.length(fdtd::AbstractFDTD) = length(fdtd.fields)
Base.ndims(fdtd::AbstractFDTD{N}) where N = N

_memorydepth(::AbstractFields{N,M}) where {N,M} = M
_memorydepth(::AbstractFDTD{N,M}) where {N,M} = M

_ndimpml(::AbstractFDTD{N,M,P}) where {N,M,P} = P
_whichdimspml(fdtd::AbstractFDTD) = _whichdimspml(fdtd.simulation)

################################################################################
# Pretty Printing

import ..Common.PRINTED_COLOR_DARK
import ..Common.PRINTED_COLOR_LIGHT
import ..Common.PRINTED_COLOR_NUMBER
import ..Common.PRINTED_COLOR_VARIABLE

function Base.show(io::IO,af::AbstractFields{N,M}) where {N,M}
    print(io,"$(N)D ")
    if typeof(af)<:HelmholtzWaveFields
        printstyled(io,"HelmholtzFields",color=PRINTED_COLOR_DARK)
    elseif typeof(af)<:HelmholtzBlochFields
        printstyled(io,"HelmholtzBlochFields",color=PRINTED_COLOR_DARK)
    end
    print(io," (")
    printstyled(io, length(af), color=PRINTED_COLOR_NUMBER)
    print(io," points, depth ")
    printstyled(io, M, color=PRINTED_COLOR_NUMBER)
    print(io,")")
end

function Base.show(io::IO,fdtd::AbstractFDTD{N,M}) where {N,M}
    print(io,"$(N)D ")
    if typeof(fdtd)<:HelmholtzFDTD
        printstyled(io,"HelmholtzFDTD",color=PRINTED_COLOR_LIGHT)
    end
    print(io," (")
    printstyled(io, length(fdtd), color=PRINTED_COLOR_NUMBER)
    print(io," points, depth ")
    printstyled(io, M, color=PRINTED_COLOR_NUMBER)
    printstyled(io,", step ")
    printstyled(io,"$(fdtd.n[1])", color=PRINTED_COLOR_NUMBER)
    print(io,", time ")
    printstyled(io,"t=",fmt("1.2e",fdtd.t[1]), color=PRINTED_COLOR_NUMBER)
    print(io,")")
end

function Base.show(io::IO, options::AnimationOptions)
    printstyled(io, "AnimationOptions: ", color=PRINTED_COLOR_DARK)
    print(io, "\n\tinterval\t")
    printstyled(io, options.interval; color=PRINTED_COLOR_NUMBER)
    print(io, "\n\tstart   \t")
    printstyled(io, options.start; color=PRINTED_COLOR_NUMBER)
    print(io, "\n\tplot by \t")
    printstyled(io, options.by; color=PRINTED_COLOR_NUMBER)
    print(io, "\n\tylims   \t")
    printstyled(IOContext(io,:type_info=>Vector), options.ylims; color=PRINTED_COLOR_NUMBER)
    print(io, "\n\toptions \t")
    printstyled(io, options.plotoptions.data, color=PRINTED_COLOR_NUMBER)
end

################################################################################
# Plotting

@recipe f(af::AbstractFields; by=abs2) = af, by
@recipe f(by::Function, af::AbstractFields) = af, by

@recipe f(sim::Simulation, af::AbstractFields; by=abs2) = sim, af, by
@recipe f(af::AbstractFields, sim::Simulation; by=abs2) = sim, af, by
@recipe f(sim::Simulation, by::Function, af::AbstractFields) = sim, af, by
@recipe f(by::Function, sim::Simulation, af::AbstractFields) = sim, af, by
@recipe f(by::Function, af::AbstractFields, sim::Simulation) = sim, af, by

@recipe f(fdtd::AbstractFDTD; by=abs2) = fdtd, by
@recipe f(by::Function, fdtd::AbstractFDTD) = fdtd, by

@recipe f(sim::Simulation, fdtd::AbstractFDTD; by=abs2) = simulation, fdtd, by
@recipe f(fdtd::AbstractFDTD, sim::Simulation; by=abs2) = simulation, fdtd, by
@recipe f(sim::Simulation, by::Function, fdtd::AbstractFDTD) = simulation, fdtd, by
@recipe f(by::Function, sim::Simulation, fdtd::AbstractFDTD) = simulation, fdtd, by
@recipe f(by::Function, fdtd::AbstractFDTD, sim::Simulation) = simulation, fdtd, by

"""
    gif(fdtd, [filename; fps])
"""
gif(fdtd::AbstractFDTD; kwargs...) = gif(fdtd.options.animation; kwargs...)

################################################################################

foreach(include,files)

end # module

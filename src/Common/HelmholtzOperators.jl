"""
    module HelmholtzOperators
"""
module HelmholtzOperators

export Helmholtz

files = (
    "1D/HelmholtzOperators1D.jl",
    "2D/Symmetric/HelmholtzOperators2D.jl",
    # "3D/HelmholtzOperators3D.jl"
    )

using ..VectorFields
using ..Simulations
using ..Dispersions
using SparseArrays

import ..Symmetric, ..Unsymmetric
import LinearAlgebra: I

"""
    maxwell(simulation{1}; ky=0, kz=0) -> m
"""
struct Helmholtz{N,CLASS,TΣ,TSIM}
    A::SparseMatrixCSC{ComplexF64,Int}
    D²::SparseMatrixCSC{ComplexF64,Int}
    αεpFχ::SparseMatrixCSC{Complex{Float64},Int}
    laplacian::SparseMatrixCSC{ComplexF64,Int}
    Σs::TΣ
    αε::SparseMatrixCSC{ComplexF64,Int}
    Fχ::SparseMatrixCSC{ComplexF64,Int}
    Fχs::Vector{SparseMatrixCSC{ComplexF64,Int}}
    dFχdψr::Array{ComplexF64,3}
    dFχdψi::Array{ComplexF64,3}
    dFχdω::Array{ComplexF64,3}
    dFχdϕ::Array{ComplexF64,3}
    simulation::TSIM
    ka::Float64
    kb::Float64
    kc::Float64
end

################################################################################
# susceptabilities and jacobians

# linear
@inline function helmholtz_susceptability!(h::Helmholtz,ω::Number)
    sim = h.sim
    @fastmath @inbounds @simd for i ∈ eachindex(h.Fχ.nzval)
        χ = susceptability(sim.χ[i], ω)
        h.Fχ.nzval[i] = sim.F[i]*χ
    end
    @fastmath @inbounds @simd for i ∈ eachindex(h.αε.nzval) h.αεpFχ.nzval[i] = h.αε.nzval[i] + h.Fχ.nzval[i] end
    return nothing
end


# nonlinear + local jacobian
@inline function helmholtz_susceptability!(h::Helmholtz, ω, ωs::Vector, ψs::ScalarField)
    sim = h.sim
    foreach(d->susceptability(d.χ,ω,ωs,ψs), sim.dispersive_domains)
    @inbounds for μ ∈ eachindex(h.Fχs)
        @inbounds for i ∈ eachindex(h.Fχ.nzval)
            if typeof(sim.χ[i])<:TwoLevelSystem
                χ::TwoLevelSystem = sim.χ[i]
                μ==1 ? h.Fχ.nzval[i] = sim.F[i]*χ.chi[i] : nothing
                h.Fχs[μ].nzval[i] = sim.F[i]*χ.chis[i,μ]
                @inbounds @simd for ν ∈ eachindex(h.Fχs)
                    h.dFχdψr[i,μ,ν] = sim.F[i]*χ.dχdψr[i,μ,ν]
                    h.dFχdψi[i,μ,ν] = sim.F[i]*χ.dχdψi[i,μ,ν]
                    h.dFχdω[i,μ,ν] = sim.F[i]*χ.dχdω[i,μ,ν]
                    h.dFχdϕ[i,μ,ν] = sim.F[i]*χ.dχdϕ[i,μ,ν]
                end
            end
        end
    end
    @fastmath @inbounds @simd for i ∈ eachindex(h.αε.nzval) h.αεpFχ.nzval[i] = h.αε.nzval[i] + h.Fχ.nzval[i] end
    return nothing
end

################################################################################
# Pretty Printing

import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_VARIABLE
import ..PRINTED_COLOR_INSTRUCTION

function Base.show(io::IO,::Helmholtz{N,CLASS}) where {N,CLASS}
    print(io,N,"D ")
    CLASS<:Symmetric ? print(io,"Symmetric ") : nothing
    CLASS<:Unsymmetric ? print(io,"Unsymmetric ") : nothing
    printstyled(io,"Helmholtz ",color=PRINTED_COLOR_DARK)
    print(io,"Operator ")
    printstyled(io,"(call w/ args ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ω",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", [",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ωs",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,",",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ψs",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,"]) -> ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"A",color=PRINTED_COLOR_VARIABLE)
end

function Base.getproperty(h::Helmholtz,sym::Symbol)
    if sym==:sim
        return getfield(h,:simulation)
    elseif sym==:ky
        return getfield(h,:kb)
    elseif sym==:kz
        return getfield(h,:kc)
    else
        return getfield(h,sym)
    end
end

function Base.propertynames(::Helmholtz{1}, private=false)
    if private
        return fieldnames(Helmholtz)
    else
        (:simulation,:ky,:kz)
    end
end

# load 1D, 2D, 3D
foreach(include,files)

end #module

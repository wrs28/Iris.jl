#NOTE: \perp does not work, reserved as infix operator, use \bigbot instead
"""
    module TwoLevelSystems

Defines two-level dispersion, and utilties for Jacobians.
"""
module TwoLevelSystems

export TwoLevelSystem

using Formatting
using LinearAlgebra
using SparseArrays
using ...ElectricFields

import ...PRINTED_COLOR_NUMBER
import ...PRINTED_COLOR_DARK
import ...INDEX_OFFSET
import ..AbstractDispersion
import ..susceptability

mutable struct TwoLevelSystem <: AbstractDispersion
    ωₐ::Float64
    D₀::Float64
    γ⟘::Float64
    h::Vector{Float64}
    hp1⁻¹::Vector{Float64}
    chi::Vector{ComplexF64}
    chis::Matrix{ComplexF64}
    dχdψr::Array{ComplexF64,3}
    dχdψi::Array{ComplexF64,3}
    dχdω::Array{ComplexF64,3}
    dχdϕ::Array{ComplexF64,3}

    TwoLevelSystem(ωₐ, D₀, γ⟘) = new(ωₐ, D₀, γ⟘)
end

function TwoLevelSystem(;ωₐ=0, D₀=0, γ⟘=1e8, kwargs...)
    ωa = get(kwargs,:ωa,ωₐ)
    omega_a = get(kwargs,:omega_a,ωa)
    omega = get(kwargs,:omega,omega_a)
    ω = get(kwargs,:ω,omega)

    γ_perp = get(kwargs,:γ_perp,γ⟘)
    γperp = get(kwargs,:γperp,γ_perp)
    γp = get(kwargs,:γp,γperp)
    gamma_perp = get(kwargs,:gamma_perp,γp)
    γ = get(kwargs,:γ,gamma_perp)

    D0 = get(kwargs,:D0,D₀)
    D = get(kwargs,:D,D0)
    return TwoLevelSystem(ω, D, γ)
end

susceptability(tls::TwoLevelSystem,ω) = χ(tls,ω)
susceptability(tls::TwoLevelSystem,args...) = χ!(tls,args...)

Base.conj(tls::TwoLevelSystem) = TwoLevelSystem(tls.ωa,tls.D0,-tls.γp)

function Base.getproperty(tls::TwoLevelSystem,sym::Symbol)
    if Base.sym_in(sym,(:γ,:γ_perp,:gamma_perp,:γp,:γperp))
        return getfield(tls,:γ⟘)
    elseif Base.sym_in(sym,(:ω,:omega_a,:ωa,:omega))
        return getfield(tls,:ωₐ)
    elseif Base.sym_in(sym,(:D,:D0))
        return getfield(tls,:D₀)
    else
        return getfield(tls,sym)
    end
end

Base.propertynames(::TwoLevelSystem) = (:γ⟘,:ωₐ,:D₀)

function Base.setproperty!(tls::TwoLevelSystem,sym::Symbol,val::Real)
    if Base.sym_in(sym,(:γ,:γ_perp,:gamma_perp,:γp,:γperp))
        return setfield!(tls,:γ⟘,float(val))
    elseif Base.sym_in(sym,(:ω,:omega_a,:ωa,:omega))
        return setfield!(tls,:ωₐ,float(val))
    elseif Base.sym_in(sym,(:D,:D0))
        return setfield!(tls,:D₀,float(val))
    else
        return setfield!(tls,sym,float(val))
    end
end

function Base.show(io::IO, tls::TwoLevelSystem)
    printstyled(io, "Two Level System",color=PRINTED_COLOR_DARK)
    print(io," (ωₐ=")
    printstyled(io, fmt("3.2f",tls.ωₐ),color=PRINTED_COLOR_NUMBER)
    print(io,", D₀=")
    printstyled(io, fmt("2.3f",tls.D₀),color=PRINTED_COLOR_NUMBER)
    print(io,", γ⟘=")
    printstyled(io, fmt("1.1f",tls.γ⟘),color=PRINTED_COLOR_NUMBER)
    print(io,")")
end

################################ AUXILLIARIES ##################################

γ(tls::TwoLevelSystem,ω) = tls.γ⟘*one(ω)/(ω-complex(tls.ωₐ,-tls.γ⟘)*one(ω))

Γ(tls::TwoLevelSystem,ω) = abs2(γ(tls,ω))

H!(tls::TwoLevelSystem) = nothing
@inline function H!(tls::TwoLevelSystem,ωs::Vector,ψs::ElectricField)
    N,M = size(ψs)
    isdefined(tls,:h) ? nothing : tls.h = zeros(Float64,N)
    isdefined(tls,:hp1⁻¹) ? nothing : tls.hp1⁻¹ = similar(tls.h)
    isdefined(tls,:chi) ? nothing : tls.chi = Vector{ComplexF64}(undef,N)
    isdefined(tls,:chis) ? nothing : tls.chis = Array{ComplexF64,2}(undef,N,M)
    isdefined(tls,:dχdψr) ? nothing : tls.dχdψr = Array{ComplexF64,3}(undef,N,M,M)
    isdefined(tls,:dχdψi) ? nothing : tls.dχdψi = Array{ComplexF64,3}(undef,N,M,M)
    isdefined(tls,:dχdω) ? nothing : tls.dχdω = Array{ComplexF64,3}(undef,N,M,M)
    isdefined(tls,:dχdϕ) ? nothing : tls.dχdϕ = Array{ComplexF64,3}(undef,N,M,M)
    size(tls.h)==(N,1) ? nothing : tls.h = Vector{Float64}(undef,N)
    size(tls.hp1⁻¹)==(N,1) ? nothing : tls.hp1⁻¹ = Vector{Float64}(undef,N)
    size(tls.chi)==(N,1) ? nothing : tls.chi = Vector{ComplexF64}(undef,N)
    size(tls.chis)==(N,M) ? nothing : tls.chis = Array{ComplexF64,2}(undef,N,M)
    size(tls.dχdψr)==(N,M) ? nothing : tls.dχdψr = Array{ComplexF64,3}(undef,N,M,M)
    size(tls.dχdψi)==(N,M) ? nothing : tls.dχdψi = Array{ComplexF64,3}(undef,N,M,M)
    size(tls.dχdω)==(N,M) ? nothing : tls.dχdω = Array{ComplexF64,3}(undef,N,M,M)
    size(tls.dχdϕ)==(N,M) ? nothing : tls.dχdϕ = Array{ComplexF64,3}(undef,N,M,M)
    ψx,ψy,ψz = ψs.x,ψs.y,ψs.z
    for i ∈ 1:N tls.h[i] = 0 end
    @inbounds for μ ∈ eachindex(ωs)
        G = Γ(tls,ωs[μ])
        @inbounds for i ∈ 1:N
            j = mod1(i,N÷3)
            tls.h[i] += G*(abs2(ψx[j,μ])+abs2(ψy[j,μ])+abs2(ψz[j,μ]))
        end
    end
    @fastmath @inbounds @simd for i ∈ 1:N tls.hp1⁻¹[i] = 1/(1+tls.h[i]) end
    return nothing
end

χ(tls::TwoLevelSystem,ω) = tls.D₀*γ(tls,ω)


@inline function χ!(tls::TwoLevelSystem,ω,ωs::Vector,ψs::ElectricField)
    index = size(ψs,1)÷2+INDEX_OFFSET
    H!(tls,ωs,ψs)
    Dg = tls.D₀*γ(tls,ω)
    gp² = tls.γp^2
    @fastmath @inbounds @simd for i ∈ eachindex(tls.chi) tls.chi[i] = Dg*tls.hp1⁻¹[i] end
    @inbounds for μ ∈ eachindex(ωs)
        Dg = tls.D₀*γ(tls,ωs[μ])
        @fastmath @inbounds @simd for i ∈ eachindex(tls.chi) tls.chis[i,μ] = Dg*tls.hp1⁻¹[i] end
    end
    ψx, ψy, ψz = ψs.x, ψs.y, ψs.z
    @inbounds for ν ∈ eachindex(ωs) # derivative wrt ν
        Gν = Γ(tls,ωs[ν])
        ων = ωs[ν]
        @inbounds for μ ∈ eachindex(ωs) # field μ
            Dgμ = tls.D₀*γ(tls,ωs[μ])
            Dgμ² = tls.D₀*γ(tls,ωs[μ])^2
            Dgμ2Gν = Dgμ*2Gν
            Dgμ2Gν² = Dgμ2Gν*Gν
            @inbounds @simd ivdep for i ∈ eachindex(tls.chi)
                k = mod1(i,length(tls.chi)÷3)
                hp1⁻² = tls.hp1⁻¹[i]^2
                tls.dχdψr[i,μ,ν], tls.dχdψi[i,μ,ν] = -Dgμ2Gν*hp1⁻².*reim(ψs[i,ν])
                tls.dχdω[i,μ,ν] = -(μ==ν)*Dgμ²*tls.hp1⁻¹[i]/tls.γp +
                    Dgμ2Gν²*hp1⁻²*(ων-tls.ωa)*(abs2(ψx[k,ν])+abs2(ψy[k,ν])+abs2(ψz[k,ν]))/gp²
                tls.dχdϕ[i,μ,ν] = -Dgμ2Gν*hp1⁻²*(abs2(ψx[k,ν])+abs2(ψy[k,ν])+abs2(ψz[k,ν]))/real(ψs[index,ν])
            end
        end
    end
    return nothing
end


"""
    mutable struct TwoLevelSystem

    TwoLevelSystem(; [ωa=0][D0=0][γperp=1e8]) -> tls
"""
TwoLevelSystem

"""
    susceptability(::TwoLevelSystem,ω) -> χ::ComplexF64
    susceptability(::TwoLevelSystem,ω,ωs,ψs)

The second computes in-place (stored in TwoLevelSystem).
"""
susceptability

end # module

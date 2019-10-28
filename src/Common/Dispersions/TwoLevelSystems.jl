#NOTE: \perp does not work, reserved as infix operator, use \bigbot instead
module TwoLevelSystems

export TwoLevelSystem

using LinearAlgebra
using SparseArrays
using ...ElectricFields

import ...PRINTED_COLOR_NUMBER
import ...PRINTED_COLOR_DARK
import ..AbstractDispersion
import ..susceptability


mutable struct TwoLevelSystem <: AbstractDispersion
    ωₐ::Float64
    D₀::Float64
    γ⟘::Float64
    h::Vector{Float64}
    hp1⁻¹::Vector{Float64}
    chi::Matrix{ComplexF64}
    dχdψr::Matrix{ComplexF64}
    dχdψi::Matrix{ComplexF64}

    TwoLevelSystem(ωₐ, D₀, γ⟘) = new(ωₐ, D₀, γ⟘)
end

function TwoLevelSystem(;ωₐ=0, D₀=0, γ⟘=1e8, kwargs...)
    ωa = get(kwargs,:ωa,ωₐ)
    omega_a = get(kwargs,:omega_a,ωa)
    ω = get(kwargs,:ω,omega_a)

    γ_perp = get(kwargs,:γ_perp,γ⟘)
    γp = get(kwargs,:γp,γ_perp)
    gamma_perp = get(kwargs,:gamma_perp,γp)
    γ = get(kwargs,:γ,gamma_perp)

    D0 = get(kwargs,:D0,D₀)
    D = get(kwargs,:D,D0)
    return TwoLevelSystem(ω, D, γ)
end

susceptability(tls::TwoLevelSystem,ω,args...) = χ(tls,ω,args...)
# susceptability(tls::TwoLevelSystem,args...) = χ(tls,args...)

Base.conj(tls::TwoLevelSystem) = TwoLevelSystem(tls.ωa,tls.D0,-tls.γp)

function Base.getproperty(tls::TwoLevelSystem,sym::Symbol)
    if Base.sym_in(sym,(:γ,:γ_perp,:gamma_perp,:γp))
        return getfield(tls,:γ⟘)
    elseif Base.sym_in(sym,(:ω,:omega_a,:ωa))
        return getfield(tls,:ωₐ)
    elseif Base.sym_in(sym,(:D,:D0))
        return getfield(tls,:D₀)
    else
        return getfield(tls,sym)
    end
end

function Base.setproperty!(tls::TwoLevelSystem,sym::Symbol,val::Real)
    if Base.sym_in(sym,(:γ,:γ_perp,:gamma_perp,:γp))
        return setfield!(tls,:γ⟘,val)
    elseif Base.sym_in(sym,(:ω,:omega_a,:ωa))
        return setfield!(tls,:ωₐ,val)
    elseif Base.sym_in(sym,(:D,:D0))
        return setfield!(tls,:D₀,val)
    else
        return setfield!(tls,sym,val)
    end
end

function Base.show(io::IO, tls::TwoLevelSystem)
    printstyled(io, "Two Level System",color=PRINTED_COLOR_DARK)
    print(io," (ωₐ=")
    printstyled(io, tls.ωₐ,color=PRINTED_COLOR_NUMBER)
    print(io,", D₀=")
    printstyled(io, tls.D₀,color=PRINTED_COLOR_NUMBER)
    print(io,", γ⟘=")
    printstyled(io, tls.γ⟘,color=PRINTED_COLOR_NUMBER)
    print(io,")")
end

################################ AUXILLIARIES ##################################

γ(tls::TwoLevelSystem,ω) = tls.γ⟘*one(ω)/(ω-complex(tls.ωₐ,-tls.γ⟘)*one(ω))

Γ(tls::TwoLevelSystem,ω) = abs2(γ(tls,ω))

@inline function H!(tls::TwoLevelSystem,ωs::Vector,ψs::ElectricField)
    N,M = size(ψs)
    isdefined(tls,:h) ? nothing : tls.h = zeros(Float64,N)
    isdefined(tls,:hp1⁻¹) ? nothing : tls.hp1⁻¹ = similar(tls.h)
    isdefined(tls,:chi) ? nothing : tls.chi = Matrix{ComplexF64}(undef,N,M)
    isdefined(tls,:dχdψr) ? nothing : tls.dχdψr = similar(chi)
    isdefined(tls,:dχdψi) ? nothing : tls.dχdψi = similar(chi)
    ψx,ψy,ψz = ψs.x,ψs.y,ψs.z
    for i ∈ 1:N tls.h[i] = 0 end
    for μ ∈ eachindex(ωs)
        G = Γ(tls,ωs[μ])
        for i ∈ 1:N
            j = mod1(i,N÷3)
            tls.h[i] += G*(abs2(ψx[j,μ])+abs2(ψy[j,μ])+abs2(ψz[j,μ]))
        end
    end
    for i ∈ 1:N tls.hp1⁻¹[i] = 1/(1+tls.h[i]) end
    return nothing
end

χ(tls::TwoLevelSystem,ω) = tls.D₀*γ(tls,ω)

@inline function χ!(tls::TwoLevelSystem,ω,ωs::Vector,ψs::ElectricField)
    H!(tls,ωs,ψs)
    Dg = tls.D₀*γ(tls,ω)
    for i ∈ eachindex(tls.h) tls.chi[i,1] = Dg*tls.hp1⁻¹[i] end
    return nothing
end

@inline function χ!(tls::TwoLevelSystem,ωs::Vector,ψs::ElectricField)
    H!(tls,ωs,ψs)
    for μ ∈ eachindex(ωs)
        ω = ωs[μ]
        Dg = tls.D₀*γ(tls,ω)
        Dg2G = Dg*2Γ(tls,ω)
        for i ∈ eachindex(tls.chi)
            tls.chi[i,μ] = Dg*tls.hp1⁻¹[i]
            tls.dχdψr[i,μ], tls.dχdψi[i,μ] = -Dg2G*(tls.hp1⁻¹[i])^2 .*reim(ψs[i,μ])
        end
    end
    return nothing
end


"""
    struct TwoLevelSystem
    TwoLevelSystem(tls; :key1 => value1, :key2 => value2, ...) -> tls
"""
TwoLevelSystem

"""
    susceptability(::TwoLevelSystem,ω)
    susceptability(::TwoLevelSystem,ωs,ψs)
"""
susceptability

end # module

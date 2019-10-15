#NOTE: \perp does not work, reserved as infix operator
# use \bigbot instead
module TwoLevelSystems

export TwoLevelSystem
# export jacobian_lasing

using LinearAlgebra
using SparseArrays

import ..AbstractDispersion
import ..susceptability


"""
    struct TwoLevelSystem
    TwoLevelSystem(tls; :key1 => value1, :key2 => value2, ...) -> tls
"""
TwoLevelSystem

"""
    susceptability(::TwoLevelSystem,ω)
    susceptability(::TwoLevelSystem,ωs,ψs)
"""

mutable struct TwoLevelSystem <: AbstractDispersion
    ωₐ::Float64
    D₀::Float64
    γ⟘::Float64
    h::Array{Float64,1}
    h⁻¹::Array{Float64,1}
    chi::Array{ComplexF64,1}

    TwoLevelSystem(ωₐ, D₀, γ⟘) = new(ωₐ, D₀, γ⟘)
end
function TwoLevelSystem(;ωₐ = 0, D₀ = 0, γ⟘ = 1e8, kwargs...)

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

Base.conj(tls::TwoLevelSystem) = TwoLevelSystem(tls.ωa,tls.D0,-tls.γp)
Base.conj(()) = ()

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
function Base.setproperty!(tls::TwoLevelSystem,sym::Symbol,val::Float64)
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
Base.setproperty!(tls::TwoLevelSystem,sym::Symbol,val::Number) = setproperty!(tls,sym,float(val))


function Base.show(io::IO, tls::TwoLevelSystem)
    print(io, "Two Level System:\n")
    print(io, "\tωₐ: ", tls.ωₐ, "\n",
    "\tD₀: ", tls.D₀, "\n",
    "\tγ⟘: ", tls.γ⟘)
end


"""
    susceptability(::TwoLevelSystem,ωs::Array,ψ::Array) -> χ::SparseMatrixCSC
"""
susceptability(tls::TwoLevelSystem,ω::Number) = χ(tls,ω)
susceptability(tls::TwoLevelSystem,args...) = χ!(tls,args...)

# function jacobian_scattering(dispersion::TwoLevelSystem,ωs::Array,ψ::Array)
#     jacobian_ψ(dispersion,ωs,ψ)
# end
#
#
# jacobian_lasing(tls::TwoLevelSystem,ωs::Array,ψ::Array) = jacobian_lasing_ψ(tls,ωs,ψ)# + jacobian_lasing_ω(tls,ωs,ψ)
#
#
# function jacobian_lasing_ψ(tls::TwoLevelSystem,ωs::Array,ψ::Array)
#     N,M = size(ψ)
#
#     x = χ(tls,ωs,ψ)
#     x_temp1 = Array{ComplexF64,1}(undef,M*N)
#     x_temp2 = Array{ComplexF64,1}(undef,M*N)
#     for μ ∈ eachindex(ωs)
#         for i ∈ 1:N
#             x_temp1[(μ-1)N+i] = -x[i,μ]*ωs[μ]^2
#             if i == N÷2 + INDEX_OFFSET
#                 x_temp2[(μ-1)N+i] = 0
#             else
#                 x_temp2[(μ-1)N+i] = -x[i,μ]*ωs[μ]^2
#             end
#         end
#     end
#     XΩ1² = spdiagm(0=>x_temp1)
#     XΩ2² = spdiagm(0=>x_temp2)
#
#     Jψr,Jψi,Jψ0 = ∂χ∂ψψ(tls,ωs,ψ)
#     J1 = Jψr+Jψ0-XΩ1²; dropzeros!(J1)
#     J2 = Jψi-XΩ2²; dropzeros!(J2)
#     Jr = hcat(real(J1),real(J2))
#     Ji = hcat(imag(J1),imag(J2))
#     return vcat(Jr,Ji)
# end
#
# function jacobian_lasing_ω(tls::TwoLevelSystem,ωs,ψ)
#     N,M = size(ψ)
#     x = χ(tls,ωs,ψ)
#     dxdwp = ∂χ∂ωψ(tls,ωs,ψ)
#
#     rows = Array{Int,1}(undef,N*M^2)
#     cols = Array{Int,1}(undef,N*M^2)
#     vals = Array{ComplexF64,1}(undef,N*M^2)
#     for μ ∈ eachindex(ωs)
#         for i ∈ 1:M*N
#             rows[(μ-1)N*M+i] = i
#             cols[(μ-1)N*M+i] = N*M+(μ-1)N+1
#             vals[(μ-1)N*M+i] = -2ωs[μ]*x[i]*ψ[i] - dxdwp[μ][i]
#         end
#     end
#     Jω = sparse(rows,cols,vals,N*M,2N*M)
#     return vcat(real(Jω),imag(Jω))
# end



##################### AUXILLIARIES FOR COMPUTING JACOBIANS #####################

χ(tls::TwoLevelSystem,ω) = tls.D₀*γ(tls,ω)

@inline function χ!(tls::TwoLevelSystem,ω,ωs::Array,ψs::Array)
    H!(tls,ωs,ψs)
    Dg = tls.D₀*γ(tls,ω)
    for i ∈ eachindex(tls.chi)
        tls.chi[i] = Dg*tls.h⁻¹[i]
    end
    return nothing
end

function χ(tls::TwoLevelSystem,ωs::Array,ψs::Array)
    N,M = size(ψs)
    chi = Array{ComplexF64,2}(undef,N,M)
    χ!(chi,tls,ωs)
    return chi
end
@inline function χ!(chi::Array{ComplexF64,2},tls::TwoLevelSystem,ωs::Array,ψs::Array)
    N,M = size(ψs)
    H!(tls,ωs,ψs)
    Dg = tls.D₀*γ(tls,ωs)
    for μ ∈ eachindex(ωs)
        for i ∈ eachindex(tls.h⁻¹)
            chi[i,μ] = Dg[μ]*tls.h⁻¹[i]
        end
    end
    return nothing
end

γ(tls::TwoLevelSystem,ωs::Array,ψs::Array) = γ(tls,ωs)
γ(tls::TwoLevelSystem,ωs::Array) = γ.(Ref(tls),ωs)
γ(tls::TwoLevelSystem,ω::Number) = tls.γ⟘/(ω-complex(tls.ωₐ,-tls.γ⟘)*one(ω))

Γ(tls::TwoLevelSystem,ωs::Array,ψs::Array) = Γ(tls,ωs)
Γ(tls::TwoLevelSystem,ωs::Array) = Γ.(Ref(tls),ωs)
Γ(tls::TwoLevelSystem,ω::Number) = abs2(γ(tls,ω))

@inline function H!(tls::TwoLevelSystem,ωs::Array,ψs::Array)
    N,M = size(ψs)
    isdefined(tls,:h) ? nothing : tls.h = zeros(Float64,N)
    isdefined(tls,:h⁻¹) ? nothing : tls.h⁻¹ = Array{Float64,1}(undef,N)
    isdefined(tls,:chi) ? nothing : tls.chi = Array{ComplexF64,1}(undef,N)
    length(tls.h)==N ? nothing : tls.h = zeros(Float64,N)
    length(tls.h⁻¹)==N ? nothing : tls.h⁻¹ = Array{Float64,1}(undef,N)
    length(tls.chi)==N ? nothing : tls.chi = Array{Float64,1}(undef,N)
    G = Γ(tls,ωs,ψs)
    for i ∈ 1:N
        tls.h[i] = 0
        for μ ∈ eachindex(ωs)
            i1 = 0N÷3 + mod1(i,N÷3)
            i2 = 1N÷3 + mod1(i,N÷3)
            i3 = 2N÷3 + mod1(i,N÷3)
            tls.h[i] += G[μ]*(abs2(ψs[i1,μ])+abs2(ψs[i2,μ])+abs2(ψs[i3,μ]))
        end
        tls.h⁻¹[i] = 1/(1+tls.h[i])
    end
    return nothing
end

# function ∂χ∂ωψ(tls::TwoLevelSystem,ωs::Array,ψ::Array)
#     N,M = size(ψ)
#     x = Array{ComplexF64,2}(undef,N,M)
#     χ!(x,tls,ωs,ψ)
#     dxdw = Array{Array{ComplexF64,2},1}(undef,M)
#     ∂χ∂ωψ!(dxdw,tls,ωs,ψ,x)
#     return map(i->dxdw[i][:].*ψ[:],eachindex(dxdw))
# end
# @inline function ∂χ∂ωψ!(dxdw::Array,tls::TwoLevelSystem,ωs::Array,ψ::Array,x::Array)
#     N,M = size(ψ)
#     h⁻¹ = 1 ./ (1 .+ H(tls,ωs,ψ))
#     dhdw = ∂H∂ω(tls,ωs,ψ)
#     g = γ(tls,ωs,ψ)
#     dg = ∂γ∂ω(tls,ωs,ψ)
#     for ρ ∈ eachindex(ωs)
#         dxdw[ρ] = Array{ComplexF64,2}(undef,N,M)
#         for μ ∈ eachindex(ωs)
#             for k ∈ eachindex(h⁻¹)
#                 dxdw[ρ][k,μ] = -x[k,μ]*h⁻¹[k]*dhdw[k,ρ]*ψ[k,μ]
#                 ρ==μ ? dxdw[μ][k,μ] += x[k,μ]*dg[μ]*ψ[k,μ]/g[μ] : nothing
#             end
#         end
#     end
#     return nothing
# end
#
# ∂γ∂ω(tls::TwoLevelSystem,ωs::Array,ψ::Array) = ∂γ∂ω(tls,ωs)
# ∂γ∂ω(tls::TwoLevelSystem,ωs::Array) = ∂γ∂ω.(Ref(tls),ωs)
# ∂γ∂ω(tls::TwoLevelSystem,ω::Number) = -γ(tls,ω)^2/tls.γ
#
# ∂Γ∂ω(tls::TwoLevelSystem,ωs::Array,ψ::Array) = ∂Γ∂ω(tls,ωs)
# ∂Γ∂ω(tls::TwoLevelSystem,ωs::Array) = ∂Γ∂ω.(Ref(tls),ωs)
# ∂Γ∂ω(tls::TwoLevelSystem,ω::Number) = -2Γ(tls,ω)^2*(ω-tls.ωₐ)/tls.γ^2
#
# function ∂H∂ω(tls::TwoLevelSystem,ωs::Array,ψ::Array)
#     a = Array{Float64,2}(undef,size(ψ))
#     ∂H∂ω!(a,tls,ωs,ψ)
#     return a
# end
# @inline function ∂H∂ω!(a::Array{Float64,2},tls::TwoLevelSystem,ωs::Array,ψ::Array)
#     b = ∂Γ∂ω(tls,ωs,ψ)
#     for μ ∈ eachindex(ωs)
#         for i ∈ 1:size(ψ,1)
#             a[i,μ] = b[μ]*abs2(ψ[i,μ])
#         end
#     end
#     return nothing
# end
#
# function ∂χ∂ψψ(tls::TwoLevelSystem,ωs::Array,ψ::Array)
#     N,M = size(ψ)
#
#     DXDPR = Array{SparseMatrixCSC{ComplexF64,Int},1}(undef,M^2)
#     DXDPI = Array{SparseMatrixCSC{ComplexF64,Int},1}(undef,M^2)
#
#     dxdpr = Array{Array{ComplexF64,2},1}(undef,M)
#     dxdpi = Array{Array{ComplexF64,2},1}(undef,M)
#     ∂χ∂ψψ!(dxdpr,dxdpi,tls,ωs,ψ)
#
#     DXDP0R = Array{SparseMatrixCSC{ComplexF64,Int},1}(undef,M^2)
#     dxdp0r = Array{Array{ComplexF64,2},1}(undef,M)
#     ∂χ∂ψ0ψ!(dxdp0r,tls,ωs,ψ)
#
#     temp0 = Array{ComplexF64,1}(undef,N)
#     temp = spdiagm(0=>temp0)
#     rows0 = Array{Int,1}(undef,N)
#     cols0 = Array{Int,1}(undef,N)
#
#     for μ ∈ eachindex(ωs)
#         for i ∈ 1:N rows0[i] = (μ-1)N+i end
#         for ν ∈ eachindex(ωs)
#             sμν = sparse([μ],[ν],ComplexF64[1],M,M)
#
#             for i ∈ 1:N temp.nzval[i] = dxdpr[ν][i,μ] end
#             DXDPR[(ν-1)M+μ] = kron(sμν,temp)
#
#             for i ∈ 1:N temp.nzval[i] = dxdpi[ν][i,μ] end
#             DXDPI[(ν-1)M+μ] = kron(sμν,temp)
#
#             for i ∈ 1:N temp0[i] = dxdp0r[ν][i,μ] end
#             for i ∈ 1:N cols0[i] = (ν-1)N+N÷2+INDEX_OFFSET end
#             DXDP0R[(ν-1)M+μ] = sparse(rows0,cols0,temp0,M*N,M*N)
#         end
#     end
#     return dropzeros!(sum(DXDPR)), dropzeros!(sum(DXDPI)), dropzeros!(sum(DXDP0R))
# end
# @inline function ∂χ∂ψψ!(dxdpr,dxdpi,tls::TwoLevelSystem,ωs::Array,ψ::Array)
#     x = Array{ComplexF64,2}(undef,size(ψ))
#     χ!(x,tls,ωs,ψ)
#     h⁻¹ = 1 ./ (1 .+ H(tls,ωs,ψ))
#     dhr, dhi = ∂H∂ψ(tls,ωs,ψ)
#
#     for ν ∈ eachindex(ωs)
#         dxdpr[ν] = Array{ComplexF64,2}(undef,size(ψ))
#         dxdpi[ν] = Array{ComplexF64,2}(undef,size(ψ))
#         for μ ∈ eachindex(ωs)
#             for k ∈ 1:size(ψ,1)
#                 if k==size(ψ,1)÷2+INDEX_OFFSET
#                     dxdpr[ν][k,μ] = 0
#                     dxdpi[ν][k,μ] = 0
#                 else
#                     dxdpr[ν][k,μ] = -x[k,μ]*h⁻¹[k]*dhr[ν][k]*ψ[k,μ]
#                     dxdpi[ν][k,μ] = -x[k,μ]*h⁻¹[k]*dhi[ν][k]*ψ[k,μ]
#                 end
#             end
#         end
#     end
#     return nothing
# end
# @inline function ∂χ∂ψ0ψ!(dxdp,tls::TwoLevelSystem,ωs::Array,ψ::Array)
#     x = Array{ComplexF64,2}(undef,size(ψ))
#     χ!(x,tls,ωs,ψ)
#     h⁻¹ = 1 ./ (1 .+ H(tls,ωs,ψ))
#
#     G = Γ(tls,ωs,ψ)
#     k = size(ψ,1)÷2+INDEX_OFFSET
#     for ν ∈ eachindex(ωs)
#         dxdp[ν] = Array{ComplexF64,2}(undef,size(ψ))
#         for μ ∈ eachindex(ωs)
#             Gν2ψ0ν = 2G[ν]/ψ[k,ν]
#             for j ∈ 1:size(ψ,1)
#                 dxdp[ν][j,μ] = -x[j,μ]*Gν2ψ0ν*abs2(ψ[j,ν])*ψ[j,μ]*h⁻¹[j]
#             end
#         end
#     end
#     return nothing
# end
#
# function ∂H∂ψ(tls::TwoLevelSystem,ωs::Array,ψ::Array)
#     N,M = size(ψ)
#     dhdpr = Array{Array{Float64,1},1}(undef,M)
#     dhdpi = Array{Array{Float64,1},1}(undef,M)
#     ∂H∂ψ!(dhdpr,dhdpi,tls,ωs,ψ)
#     return dhdpr, dhdpi
# end
# @inline function ∂H∂ψ!(dhdpr::Array,dhdpi::Array,tls::TwoLevelSystem,ωs::Array,ψ::Array)
#     G = Γ(tls,ωs,ψ)
#     for ρ ∈ eachindex(ωs)
#         dhdpr[ρ] = Array{Float64,1}(undef,size(ψ,1))
#         dhdpi[ρ] = Array{Float64,1}(undef,size(ψ,1))
#         for k ∈ 1:size(ψ,1)
#             dhdpr[ρ][k] = 2G[ρ]*real(ψ[k,ρ])
#             dhdpi[ρ][k] = 2G[ρ]*real(ψ[k,ρ])
#         end
#     end
#     return nothing
# end


end # module

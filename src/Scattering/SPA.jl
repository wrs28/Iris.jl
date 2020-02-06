"""
*S*ingle *P*ole *A*pproximation for extra nonlinear scattering solutions
"""
module SPA

export spa

using ...Common
using LinearAlgebra
using NLsolve
using ProgressMeter
using RecipesBase

import ..HelmholtzNLS
import ..HelmholtzCF
import ..helmholtzeigen
import ..scattering!
import ...Common.SPA_ACCEPTANCE_PHASE

################################################################################
# spa main

"""
    spa(nls; [refine=false, nev=1, verbose=false])

Single Pole Approximation of nonlinear scattering problem `nls` for solutions besides
that given in `nls.ψ`.

`refine` | if `true`, additional refinement through nonlinear CF equation
`nev=1` | number of CF states used in refinement, making it a "few pole approximation".
`verbose` | if `true`, shows traces for CF nonlinear refinement
"""
function spa(nls::HelmholtzNLS; refine::Bool=false, forcerefine::Bool=false, nev::Integer=1, kwargs...)
    sim = nls.simulation

    # generate CF states
        cf = HelmholtzCF(sim)
        η, u = helmholtzeigen(cf, nls.ω, [nls.ω], nls.ψ; nev=nev)

    # initialize CFAmplitude structure with CF eigenpairs
        cfa1 = CFAmplitudeVector(nls, η, u)
        cfa2 = CFAmplitudeVector(nls, η, u)
        cfa3 = CFAmplitudeVector(nls, η, u)
        cfa4 = CFAmplitudeVector(nls, η, u)

    # check that nonlinear solution exists
    if nls.converged[]
        # single-pole approximations
            phase = (angle(sum(sim.F.* abs2.(u[:,1]) .* u[:,1] .* nls.ψ[:])) - angle(sum(sim.F.*u[:,1].^3 .*conj(nls.ψ[:]))))/2
            τ = nls.ψ[:]*cis(-phase)

            A = η[1]*cfa1.Γ*sum(sim.F .* abs2.(u[:,1]) .* u[:,1].^2)*sim.dx
            B = (η[1] + cfa1.γ*cfa1.D₀)*cfa1.Γ*sum(sim.F .* abs2.(u[:,1]) .*u[:,1] .* τ)*sim.dx
            C = η[1]*cfa1.Γ*sum(sim.F .* u[:,1].^3 .*τ)*sim.dx
            E = cfa1.γ*cfa1.D₀*cfa1.Γ*sum(sim.F .* abs2.(u[:,1]) .* τ.^2)*sim.dx
            G = (η[1] - cfa1.γ*cfa1.D₀) + η[1]*cfa1.Γ*sum(sim.F .* u[:,1].^2 .* τ.^2)*sim.dx

            Σ = (B+C)/A
            Π = (E+G)/A
            cfa1.Σ, cfa1.Π = Σ, Π
            cfa2.Σ, cfa2.Π = Σ, Π

        small_phaseΣ = abs(angle(Σ)) < SPA_ACCEPTANCE_PHASE || abs(angle(Σ)-π) < SPA_ACCEPTANCE_PHASE || abs(angle(Σ)+π) < SPA_ACCEPTANCE_PHASE
        small_phaseΠ = abs(angle(Π)) < SPA_ACCEPTANCE_PHASE || abs(angle(Π)-π) < SPA_ACCEPTANCE_PHASE || abs(angle(Π)+π) < SPA_ACCEPTANCE_PHASE
        if small_phaseΣ && small_phaseΠ
            cfa1.multivalued[] = real(Σ/2)^2 > real(Π)
            cfa2.multivalued[] = real(Σ/2)^2 > real(Π)
        else
            cfa1.multivalued[] = false
            cfa2.multivalued[] = false
        end
        cfa1.estimate = -cis(phase)*abs( (Σ/2) + sqrt((Σ/2)^2 - Π) )
        cfa2.estimate = -cis(phase)*abs( (Σ/2) - sqrt((Σ/2)^2 - Π) )

        if forcerefine || (refine && cfa1.multivalued[])
            #cfa1
                cfa1.refined[] = true
                cfa2.refined[] = true
                x = a_to_x(cfa1,nev)
                res = nlsolve(cfa1, x; kwargs...)
                cfa1.converged[] = converged(res)
                x_to_a!(cfa1,res.zero)

            append!(cfa2.a0,cfa1.a)

            #cfa2
                x = a_to_x(cfa2,nev)
                res = nlsolve(cfa2, x; kwargs...)
                cfa2.converged[] = converged(res)
                x_to_a!(cfa2,res.zero)
        end

        ρ = sqrt((E-G)/(A+C*(C-B)/E))
        θ = phase + acos((C-B)*ρ/2/E)
        cfa3.Σ, cfa3.Π = ρ, θ
        θ = phase - acos((C-B)*ρ/2/E)
        cfa4.Σ, cfa4.Π = ρ, θ

        small_phaseρ = abs(angle(ρ)) < SPA_ACCEPTANCE_PHASE || abs(angle(ρ)-π) < SPA_ACCEPTANCE_PHASE || abs(angle(ρ)+π) < SPA_ACCEPTANCE_PHASE
        small_phaseθ = abs((C-B)*ρ/2/E) ≤ 1#abs(angle(θ)) < SPA_ACCEPTANCE_PHASE || abs(angle(θ)-π) < SPA_ACCEPTANCE_PHASE || abs(angle(θ)+π) < SPA_ACCEPTANCE_PHASE
        if small_phaseρ && small_phaseθ
            cfa3.multivalued[] = true
            cfa4.multivalued[] = true
        else
            cfa3.multivalued[] = false
            cfa4.multivalued[] = false
        end
        cfa3.estimate = cis(real(cfa3.θ))*real(cfa3.ρ)
        cfa4.estimate = cis(real(cfa4.θ))*real(cfa4.ρ)

        if forcerefine || (refine && cfa3.multivalued[])
            #cfa3
                cfa3.refined[] = true
                cfa4.refined[] = true

                x = a_to_x(cfa3,nev)
                res = nlsolve(cfa3, x; kwargs...)
                cfa3.converged[] = converged(res)
                x_to_a!(cfa3,res.zero)

            append!(cfa4.a0,cfa3.a)

            #cfa4
                x = a_to_x(cfa4,nev)
                res = nlsolve(cfa4, x; kwargs...)
                cfa4.converged[] = converged(res)
                x_to_a!(cfa4,res.zero)
        end
    elseif get(kwargs,:show_trace,false)
        @warn "must solve NonLinearScatteringProblem first"
    end
    return cfa1, cfa2, cfa3, cfa4
end

################################################################################
# structures

# structure used in nonlinear CF solve
struct CFAmplitudeVector{TF}
    estimate::Base.RefValue{ComplexF64}
    Σ::Base.RefValue{ComplexF64}
    Π::Base.RefValue{ComplexF64}
    a::Vector{ComplexF64}
    a0::Vector{ComplexF64}
    φ::TF
    ζ::TF
    u::TF
    η::Vector{ComplexF64}
    B::Vector{ComplexF64}
    dx::Float64
    D₀::Float64
    γ::ComplexF64
    Γ::Float64
    F::Vector{Float64}
    converged::Base.RefValue{Bool}
    multivalued::Base.RefValue{Bool}
    refined::Base.RefValue{Bool}
end

"""
    CFAmplitudeVector(nls, η, u, [a0]) -> cfa

construct nonlinear CF solver object with CF evals `η` and efuns `u`.
`estimate` is the SPA estimate, `a0` is a previously found solution to be avoided
"""
function CFAmplitudeVector(
            nls::HelmholtzNLS{N},
            η::Vector,
            u::ScalarField{N},
            a0::Vector=ComplexF64[]
            ) where N

    nev = size(u,2)
    sim = nls.simulation
    φ = nls.solution.total

    tls = nls.simulation.dispersive_domains[1].χ
    Γ = Common.Dispersions.TwoLevelSystems.Γ(tls,nls.ω)
    γ = Common.Dispersions.TwoLevelSystems.γ(tls,nls.ω)
    D₀ = tls.D₀

    a = zeros(ComplexF64,nev)

    ζ = ScalarField(sim,1)
    B = Vector{ComplexF64}(undef, nev)
    Fφ = sim.F .*φ[:] ./ (1 .+ Γ .* abs2.(φ[:]) )
    @inbounds for i ∈ eachindex(B) B[i] = γ*D₀*sum(Fφ .* u[:,i]) * sim.dx end
    return CFAmplitudeVector(Ref(complex(0.0)), Ref(complex(0.0)), Ref(complex(0.0)), a, a0, φ, ζ, u, η, B, sim.dx, D₀, γ, Γ, sim.F, Ref(true), Ref(false), Ref(false))
end

# create wavefunction φ + ∑ⱼaⱼuⱼ
"""
    (cfa)([nev]) -> φ

Generate `φ` from background field in `cfa` plus CF corrections
"""
function (cfa::CFAmplitudeVector)(nev::Integer=length(cfa.a))
    copyto!(cfa.ζ.values,cfa.φ.values)
    @inbounds for i ∈ 1:nev
        @inbounds @fastmath @simd for j ∈ eachindex(cfa.ζ)
            cfa.ζ[j] = cfa.ζ[j] + cfa.a[i]*cfa.u[j,i]
        end
    end
    return cfa.ζ
end

# cfa(F,x) used in nlsolve.
function (cfa::CFAmplitudeVector)(F::Vector{Float64},x::Vector{Float64})
    n = length(x)÷2
    # load data
        x_to_a!(cfa,x)
    # create wavefunction associated with a's
        ψ = cfa(n)
    # normalization penalizes being close to 0 or a0
        𝒩⁻¹ = 1/norm(cfa.a) + (isempty(cfa.a0) ? 0.0 : 1/norm(cfa.a-cfa.a0) )
    @inbounds @fastmath for i ∈ 1:n
        A = complex(0.0)
        @inbounds @fastmath @simd ivdep for j ∈ eachindex(cfa.F)
            A += cfa.F[j] * ψ.values[j] * cfa.u.values[j,i] / (1 + cfa.Γ * abs2(ψ.values[j]) )
        end
        A *= cfa.γ*cfa.D₀*cfa.dx
        F[i], F[n+i] = reim((-cfa.η[i]*cfa.a[i] + A - cfa.B[i])*𝒩⁻¹)
    end
    return nothing
end


function Base.getproperty(cfa::CFAmplitudeVector, sym::Symbol)
    if Base.sym_in(sym,(:Σ,:σ,:sigma,:Sigma,:s,:S,:ρ,:rho,:r))
        return getfield(cfa,:Σ)[]
    elseif Base.sym_in(sym,(:Π,:π,:Pi,:pi,:P,:p,:θ,:ϑ,:theta))
        return getfield(cfa,:Π)[]
    elseif sym == :estimate
        return getfield(cfa,:estimate)[]
    else
        return getfield(cfa,sym)
    end
end

function Base.setproperty!(cfa::CFAmplitudeVector, sym::Symbol, x::Number)
    if sym == :estimate
        getfield(cfa,:estimate)[] = x
        return cfa.a[1] = x
    elseif Base.sym_in(sym, (:Σ,:Π))
        getfield(cfa,sym)[] = x
    else
        return setfield!(cfa,sym,x)
    end
end


################################################################################
# data loaders

a_to_x(cfa::CFAmplitudeVector,n::Integer) = a_to_x(cfa.a,n)
a_to_x(a::Vector,n::Integer) = a_to_x!(Vector{Float64}(undef,2n),a)

a_to_x!(x::Vector, cfa::CFAmplitudeVector) = a_to_x!(x,cfa.a)
@inline function a_to_x!(x::Vector,a::Vector)
    n = length(x)÷2
    @inbounds @simd ivdep for i ∈ 1:n x[i], x[n+i] = reim(a[i]) end
    return x
end

x_to_a!(cfa::CFAmplitudeVector, x::Vector) = x_to_a!(cfa.a,x)
@inline function x_to_a!(a::Vector,x::Vector)
    n = length(x)÷2
    @inbounds @simd ivdep for i ∈ 1:n a[i] = complex(x[i], x[n+i]) end
    return a
end


################################################################################
# Pretty Printing


import ..Common.PRINTED_COLOR_LIGHT
import ..Common.PRINTED_COLOR_GOOD
import ..Common.PRINTED_COLOR_NUMBER
import ..Common.PRINTED_COLOR_BAD
import ..Common.PRINTED_COLOR_WARN
import ..Common.PRINTED_COLOR_VARIABLE
import ..Common.PRINTED_COLOR_INSTRUCTION

function Base.show(io::IO,cfa::CFAmplitudeVector)
    printstyled(io,"CF Amplitudes ",color=PRINTED_COLOR_LIGHT)
    if cfa.refined[]
        if cfa.converged[]
            printstyled(io,"converged ", color=PRINTED_COLOR_GOOD)
        else
            printstyled(io,"not converged ", color=PRINTED_COLOR_BAD)
        end
    end
    if cfa.multivalued[]
        printstyled(io,"multi-valued ", color=PRINTED_COLOR_WARN)
    else
        printstyled(io,"single-valued ", color=PRINTED_COLOR_INSTRUCTION)
    end
    show(IOContext(io,:compact=>true),cfa.a)
end

function Base.show(io::IO,mime::MIME"text/plain",cfa::CFAmplitudeVector)
    printstyled(io,"CF Amplitudes ",color=PRINTED_COLOR_LIGHT)
    if cfa.refined[]
        if cfa.converged[]
            printstyled(io,"converged ", color=PRINTED_COLOR_GOOD)
        else
            printstyled(io,"not converged ", color=PRINTED_COLOR_BAD)
        end
    end
    if cfa.multivalued[]
        printstyled(io,"multi-valued ", color=PRINTED_COLOR_WARN)
    else
        printstyled(io,"single-valued ", color=PRINTED_COLOR_INSTRUCTION)
    end
    print(io,"(SPA estimate of ")
    printstyled(io,"a[1]",color=PRINTED_COLOR_VARIABLE)
    print(io,": ")
    printstyled(IOContext(io,:compact=>true), cfa.estimate, color=PRINTED_COLOR_NUMBER)
    print(io,")")
    print(io,"\nVector ")
    printstyled(io,"a",color=PRINTED_COLOR_VARIABLE)
    print(io,": ")
    show(IOContext(io,:compact=>true),mime,cfa.a)
end


################################################################################
# Plotting

@recipe function f(cfa::CFAmplitudeVector)
    cfa()
end

end # module


using .SPA

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
function spa(nls::HelmholtzNLS; refine::Bool=false, forcerefine::Bool=false, nev::Integer=1, takemeout=false, nev_th::Integer=100, kwargs...)
    sim = nls.simulation

    # generate CF states
        cf = HelmholtzCF(sim)
        if nev > nev_th
            u = ScalarField(sim, 2nev_th + nev + 1)
            N = nev_th
            η1, u1 = helmholtzeigen(cf, nls.ω, [nls.ω], nls.ψ; nev=nev_th)
            η = fill(η1[1],size(u,2))
            println(length(η))
            for i ∈ 1:N u.values[:,i] = u1[:,i] end
            for i ∈ 1:N η[i] = η1[i] end
            F = Vector(diag(cf.F))
            while nev > N
                η1, u1 = helmholtzeigen(cf, nls.ω, [nls.ω], nls.ψ; η=minimum(real(η)), nev=nev_th)
                η2, u2 = helmholtzeigen(cf, nls.ω, [nls.ω], nls.ψ; η=maximum(real(η)), nev=nev_th)
                dinds1 = Int[]
                dinds2 = Int[]
                for i ∈ 1:N
                    for j ∈ eachindex(η1)
                        if abs(η[i] - η1[j]) < 1e-1
                            if abs(1-abs(sum(u.values[:,i].*u1.values[:,j].*F)*sim.dx)) < 1e-1
                                push!(dinds1,j)
                            end
                        end
                    end
                    for j ∈ eachindex(η2)
                        if abs(η[i] - η2[j]) < 1e-1
                            if abs(1-abs(sum(u.values[:,i].*u2.values[:,j].*F)*sim.dx)) < 1e-1
                                push!(dinds2,j)
                            end
                        end
                    end
                end
                kinds1 = setdiff(1:nev_th,dinds1)
                kinds2 = setdiff(1:nev_th,dinds2)
                sort!(dinds1)
                sort!(dinds2)
                deleteat!(η1,dinds1)
                deleteat!(η2,dinds2)
                N1 = length(η1)
                N2 = length(η2)
                @show size(η)
                @show size(N)
                @show size(N1)
                @show size(N2)
                for i ∈ 1:N1 η[N + i] = η1[i] end
                for i ∈ 1:N2 η[N + N1 + i] = η2[i] end
                for i ∈ 1:N1 u.values[:,N+i] = u1.values[:,kinds1[i]] end
                for i ∈ 1:N2 u.values[:,N+N1+i] = u2.values[:,kinds2[i]] end
                N += N1 + N2
            end
            η = η[1:N]
            perm = sortperm(η,by=abs)
            u0 = u
            u = ScalarField(sim,nev)
            for i ∈ 1:nev u.values[:,i] = u0.values[:,perm[i]] end
            η = η[perm[1:nev]]
        else
            η, u = helmholtzeigen(cf, nls.ω, [nls.ω], nls.ψ; nev=nev)
            F = Vector(diag(cf.F))
        end

    # initialize CFAmplitude structure with CF eigenpairs
        cfa1 = CFAmplitudeVector(nls, η, u)
        cfa2 = CFAmplitudeVector(nls, η, u)
        cfa3 = CFAmplitudeVector(nls, η, u)
        cfa4 = CFAmplitudeVector(nls, η, u)

    # check that nonlinear solution exists
    if nls.converged[]
        # single-cf approximations
            phase = angle(sum(F.* abs2.(u[:,1]) .* u[:,1] .* nls.solution.total[:]))
            τ = cis(-phase)*nls.solution.total[:]

            A = η[1]*cfa1.Γ*sum(F .* abs2.(u[:,1]) .* u[:,1].^2)*sim.dx
            B = (η[1] + cfa1.γ*cfa1.D₀)*cfa1.Γ*sum(F .* abs2.(u[:,1]) .* u[:,1] .* τ)*sim.dx
            C = η[1]*cfa1.Γ*sum(F .* u[:,1].^3 .* conj(τ))*sim.dx
            E = cfa1.γ*cfa1.D₀*cfa1.Γ*sum(F .* abs2.(u[:,1]) .* τ.^2)*sim.dx
            G = -cfa1.γ*cfa1.D₀ + η[1]*sum(sim.F .* u[:,1].^2)*sim.dx # yes, this one supposed to be sim.F

            small_angle_condition = abs(sin(angle(G/A))) ≤ abs(sin(SPA_ACCEPTANCE_PHASE))

            b = real(B/A)
            c = real(C/A)
            e = real(E/A)
            g = real(G/A)

            Σ = (b+c); cfa1.Σ = Σ; cfa2.Σ = Σ
            Π = (e+g); cfa1.Π = Π; cfa2.Π = Π
            Δ = (Σ/2)^2 - Π; cfa1.Δ = Δ; cfa2.Δ = Δ

            cfa1.scf_estimate = -cis(phase)*(Σ/2 - sqrt(complex(Δ)))
            cfa2.scf_estimate = -cis(phase)*(Σ/2 + sqrt(complex(Δ)))

            # can't exect to resolve better than imaginary parts of these, since in theory they should vanish
            discriminant_condition = Δ > 2maximum(abs∘imag,[(B+C)/A,(E+G)/A])

        if discriminant_condition && small_angle_condition
            cfa1.multivalued = true
            cfa2.multivalued = true
        else
            cfa1.multivalued = false
            cfa2.multivalued = false
        end

        # single pole approximations
            # cfa1
            ρ, θ = abs(cfa1.scf_estimate), angle(cfa1.scf_estimate)
            spa1 = SPAsolver(cfa1)
            try # try spa, if works record results, otherwise declare single-valued
                res = nlsolve(spa1,[ρ,θ]; method=:anderson, m=1, iterations=100, ftol=1e-6)
                if converged(res)
                    cfa1.converged = true
                    cfa1.multivalued = true
                    cfa1.spa_estimate = res.zero[1]*cis(res.zero[2])
                else
                    cfa1.converged = false
                    cfa1.multivalued = false
                end
            catch
                cfa1.converged = false
                cfa1.multivalued = false
            end
            if takemeout
            function f!(Fg,x)
                a = Vector{ComplexF64}(undef,length(x)÷2)
                for i ∈ eachindex(a) a[i] = complex(x[2(i-1)+1],x[2(i-1)+2]) end
                fill!(cfa1.a,0)
                copyto!(cfa1.a,a)
                ψ = cfa1(length(a))
                # Fg[1],Fg[2] = reim(A*abs2(a)+cis(phase)*B*conj(a)+cis(-phase)*C*a+cis(2phase)*E*conj(a)/a+G)
                hψ = 1 .+ cfa1.Γ .* abs2.(ψ.values)
                hφ = 1 .+ cfa1.Γ .* abs2.(cfa1.φ.values)
                q = zeros(ComplexF64,size(u,1))
                for i ∈ eachindex(a) q += a[i]*η[i]*cfa1.u[:,i] end
                for i ∈ eachindex(a)
                    A = sum(sim.F .* ψ.values  .* cfa1.u.values[:,i])*sim.dx*cfa1.γ*cfa1.D₀
                    B = sum(F .* cfa1.φ[:] .* cfa1.u.values[:,i] .* hψ)*sim.dx*cfa1.γ*cfa1.D₀
                    H = sum(hψ .* F.* q .* u[:,i])*sim.dx
                    Fg[2(i-1)+1], Fg[2(i-1)+2] = reim((-H+A-B)/norm(a))
                end
            end
            x = Vector{Float64}(undef,2nev)
            for i ∈ 1:nev x[2(i-1)+1], x[2(i-1)+2] = reim(cfa1.a[i]) end
            return nlsolve(f!,x)
        end

            #cfa2
            ρ, θ = abs(cfa2.scf_estimate), angle(cfa2.scf_estimate)
            spa2 = SPAsolver(cfa2)
            try # try spa, if works record results, otherwise declare single-valued
                res = nlsolve(spa2,[ρ,θ]; method=:anderson, m=1, iterations=100, ftol=1e-6)
                if converged(res)
                    cfa2.converged = true
                    cfa2.multivalued = true
                    cfa2.spa_estimate = res.zero[1]*cis(res.zero[2])
                else
                    cfa2.converged = false
                    cfa2.multivalued = false
                end
            catch
                cfa2.converged = false
                cfa2.multivalued = false
            end

        if forcerefine || (refine && cfa1.multivalued)
            #cfa1
                cfa1.refined = true
                x = a_to_x(cfa1,nev)
                res = nlsolve(cfa1, x; kwargs...)
                cfa1.converged = converged(res)
                x_to_a!(cfa1,res.zero)

            converged(res) ? append!(cfa2.a0,cfa1.a) : nothing

            #cfa2
                cfa2.refined = true
                x = a_to_x(cfa2,nev)
                res = nlsolve(cfa2, x; kwargs...)
                cfa2.converged = converged(res)
                x_to_a!(cfa2,res.zero)
        end

        ρ² = (e-g)/(1+c*(c-b)/e)
        if ρ² > 0
            ρ = sqrt(ρ²)
            cosθ = (c-b)*ρ/2/e
            if cosθ^2 ≤ 1
                θ = phase + acos(cosθ)
                cfa3.Σ, cfa3.Π = ρ, θ
                θ = phase - acos(cosθ)
                cfa4.Σ, cfa4.Π = ρ, θ
                cfa3.scf_estimate = cis(real(cfa3.θ))*real(cfa3.ρ)
                cfa4.scf_estimate = cis(real(cfa4.θ))*real(cfa4.ρ)
                cfa3.multivalued = true
                cfa4.multivalued = true
            else
                cfa3.multivalued = false
                cfa4.multivalued = false
            end
        else
            cfa3.multivalued = false
            cfa4.multivalued = false
        end

        if forcerefine || (refine && cfa3.multivalued)
            #cfa3
                cfa3.refined = true
                x = a_to_x(cfa3,nev)
                res = nlsolve(cfa3, x; kwargs...)
                cfa3.converged = converged(res)
                x_to_a!(cfa3,res.zero)

            converged(res) ? append!(cfa4.a0,cfa3.a) : nothing

            #cfa4
                cfa4.refined = true
                x = a_to_x(cfa4,nev)
                res = nlsolve(cfa4, x; kwargs...)
                cfa4.converged = converged(res)
                x_to_a!(cfa4,res.zero)
        end
    else #if get(kwargs,:show_trace,false)
        @warn "must solve NonLinearScatteringProblem first"
    end

    return cfa1, cfa2, cfa3, cfa4
end


################################################################################
# structures

# structure used in nonlinear CF solve
struct CFAmplitudeVector{TF}
    scf_estimate::Base.RefValue{ComplexF64}
    spa_estimate::Base.RefValue{ComplexF64}
    Σ::Base.RefValue{Float64}
    Π::Base.RefValue{Float64}
    Δ::Base.RefValue{Float64}
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
    Fφ = sim.F .* φ[:] ./ (1 .+ Γ .* abs2.(φ[:]) )
    @inbounds for i ∈ eachindex(B) B[i] = γ*D₀*sum(Fφ .* u[:,i]) * sim.dx end
    return CFAmplitudeVector(Ref(complex(0.0)), Ref(complex(0.0)), Ref(0.0), Ref(0.0), Ref(0.0), a, a0, φ, ζ, u, η, B, sim.dx, D₀, γ, Γ, sim.F, Ref(true), Ref(false), Ref(false))
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
    elseif Base.sym_in(sym,(:Δ,:δ,:d,:D,:delta,:Delta,:discriminant,:Discriminant,:dis,:Dis))
        return getfield(cfa,:Δ)[]
    elseif Base.sym_in(sym, (:scf_estimate, :spa_estimate, :converged, :multivalued, :refined))
        return getfield(cfa,sym)[]
    else
        return getfield(cfa,sym)
    end
end

function Base.setproperty!(cfa::CFAmplitudeVector, sym::Symbol, x::Number)
    if Base.sym_in(sym,(:scf_estimate, :spa_estimate))
        getfield(cfa,sym)[] = x
        return cfa.a[1] = x
    elseif Base.sym_in(sym, (:Σ, :Π, :Δ, :converged, :multivalued, :refined))
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

struct SPAsolver
    Γ::Float64
    φ::Vector{ComplexF64}
    u::Vector{ComplexF64}
    η::ComplexF64
    F::Vector{Float64}
    dx::Float64
    γ::ComplexF64
    D₀::Float64
    h::Vector{Float64}
end

function SPAsolver(cfa::CFAmplitudeVector)
    F = cfa.F ./ (1 .+ cfa.Γ .* abs2.(cfa.φ[:]) )
    SPAsolver(cfa.Γ, cfa.φ.values[:], cfa.u.values[:,1], cfa.η[1], F, cfa.dx, cfa.γ, cfa.D₀, Vector{Float64}(undef,length(F)))
end

function (spa::SPAsolver)(F::Vector,x::Vector)
    a = x[1]*cis(x[2])
    for i ∈ eachindex(spa.h) spa.h[i] = 1/(1 + spa.Γ*abs2( spa.φ[i] + a*spa.u[i])) end
    I1 = sum(spa.F .* spa.u.^2 .* spa.h)*spa.dx
    I2 = sum(spa.F .* abs2.(spa.u) .* spa.φ.^2 .* spa.h)*spa.dx*spa.Γ
    I3 = sum(spa.F .* abs2.(spa.u) .* spa.u .* spa.φ .* spa.h)*spa.dx*spa.Γ
    b, β = reim((I1 - spa.η/(spa.γ*spa.D₀))/I3)
    c, κ = reim(-I2/I3)
    θ = -atan(β+κ,b-c)
    ρ = (b+c)*cos(θ) - (β-κ)*sin(θ)
    F[1], F[2] = ρ-x[1], θ-x[2]
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
        show(IOContext(io,:compact=>true),cfa.a)
    else
        printstyled(io,"single-valued", color=PRINTED_COLOR_INSTRUCTION)
    end
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
        print(io,"(SPA estimate of ")
        printstyled(io,"a[1]",color=PRINTED_COLOR_VARIABLE)
        print(io,": ")
        printstyled(IOContext(io,:compact=>true), cfa.spa_estimate, color=PRINTED_COLOR_NUMBER)
        print(io,")")
        print(io,"\nVector ")
        printstyled(io,"a",color=PRINTED_COLOR_VARIABLE)
        print(io,": ")
        show(IOContext(io,:compact=>true),mime,cfa.a)
    else
        printstyled(io,"single-valued", color=PRINTED_COLOR_INSTRUCTION)
    end
end


################################################################################
# Plotting

@recipe function f(cfa::CFAmplitudeVector)
    cfa()
end

end # module


using .SPA

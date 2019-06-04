"""
	module Spectral
"""
module Spectral

export resonance_eigenproblem,
cf_eigenproblem,
resonance_nonlinear_eigenproblem,
eig_kl,
eig_cf,
eig_knl

using ..Defaults
using ..IrosBase

using ArnoldiMethod
using ArnoldiMethodTransformations
using Distributed
using LinearAlgebra
using NonlinearEigenproblems
using ProgressMeter
using Random
using RecipesBase
using SparseArrays
using Statistics


include("linear.jl")
include("nonlinear.jl")
include("contour_beyn_progress_bar.jl")


"""
    eig_kl(sim, k, [ka=0, kb=0; verbose, lupack, kwargs...]) -> k, ψ

Linear eigenvalue solve for `k`.

Keyword `verbose` defaults to `false`.
Keyword `lupack` defaults to `:auto` and contrls which package is used in the eigensolver.
See docs for ArnoldiMethodWrapper for details.
"""
function eig_kl(sim::Simulation, k::Number, ka::Number=0, kb::Number=0; verbose::Bool=false, lupack::Symbol=:auto, kwargs...)
    A, B, σ = resonance_eigenproblem(sim, k, ka, kb)
    λ, v, history = partialeigen(A, B, σ; diag_inv_B=true, lupack=lupack, kwargs...)

    @assert history.converged history
    verbose ? println(history) : nothing

    λ = sqrt.(λ)
    # Normalize wavefunctions according to (ψ₁,ψ₂)=δ₁₂, which requires transformed ε or F
    normalize!(sim,v,B)
    return λ, v
end


"""
    eig_cf(sim, k, [ka=0, kb=0; η=0, verbose, lupack, kwargs...]) -> η, u

Linear CF eigenvalues closest to `η`

Keyword `verbose` defaults to `false`.
Keyword `lupack` defaults to `:auto` and contrls which package is used in the eigensolver.
See docs for ArnoldiMethodWrapper for details.
"""
function eig_cf(sim::Simulation, k::Number, ka::Number=0, kb::Number=0; η::Number=0, verbose::Bool=false, lupack::Symbol=:auto, kwargs...)
    A, B, σ = cf_eigenproblem(sim, k, ka, kb)

    λ, v, history = partialeigen(A, B, η; diag_inv_B=true, lupack=lupack, kwargs...)
    @assert history.converged history
    verbose ? println(history) : nothing

    # Normalize wavefunctions according to (ψ₁,ψ₂)=δ₁₂, which requires transformed ε or F
    normalize!(sim,v,B)
    return λ, v
end


"""
    eig_knl(sim, k, ka=0, kb=0; method=contour_beyn, nk=3, display=false, quad_n=100, kwargs...) -> k,ψ
"""
function eig_knl(sim::Simulation, k::Number, ka::Number=0, kb::Number=0;
        quad_n::Int=100,
        display::Bool=false,
        method::Function=contour_beyn,
        nev::Int=3,
        quad_method=nprocs()>1 ? :ptrapz_parallel : :ptrapz,
        kwargs...
        )

    nep = resonance_nonlinear_eigenproblem(sim, k, ka, kb; check_consistency=false)
    displaylevel = display ? 1 : 0
    if display && method==contour_beyn
        k, ψ = contour_beyn(nep, true; N=quad_n, σ=k, quad_method=quad_method, neigs=nev, kwargs...)
    else
        k, ψ = method(nep; N=quad_n, σ=k, displaylevel=displaylevel, neigs=nev, quad_method=quad_method, kwargs...)
    end
    return k, ψ
end


function LinearAlgebra.normalize!(sim::Simulation,ψ,B)
    # dx = sim.dis.dx
    for i ∈ 1:size(ψ,2)
        𝒩² = sum((ψ[:,i].^2).*diag(B))#*(isinf(dx[1]) ? 1 : dx[1])*(isinf(dx[2]) ? 1 : dx[2])
        ψ[:,i] /= sqrt(𝒩²)*exp(complex(0,angle(ψ[end÷2-1,i])))
    end
    return nothing
end


@recipe function f(sim::Simulation,ψ::Array,inds=1:size(ψ,2); by=nothing, structure=false)
	@assert issubset(inds,1:size(ψ,2)) "indices $inds inconsistent with size(ψ,2)=$(size(ψ,2))"
	legend --> false
	aspect_ratio --> 1
	n = length(inds)
    if isnothing(by)
        if structure
            layout --> (1+n,3)
			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/sqrt(9+n^2+1)
        else
            layout --> (n,3)
			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/sqrt(9+n^2)
        end
    else
        if structure
            layout --> (3+n)
			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/sqrt(9+n^2)
        else
            layout --> (n)
			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/n
        end
    end
	if structure
		if isnothing(by)
			@series (sim,by)
		else
			@series (sim,real)
		end
	end
	for i ∈ 1:n
		if isnothing(by)
			bys = [:real, :imag, :abs2]
			for j ∈ eachindex(bys)
				if structure
					@series begin
						subplot --> 3i+j
						(sim,ψ,inds[i],bys[j])
					end
				else
					@series begin
						subplot --> 3(i-1)+j
						(sim,ψ,inds[i],bys[j])
					end
				end
			end
		else
			@series begin
				subplot --> i
				colorbar --> false
				(sim,ψ,inds[i],by)
			end
		end
	end
end
@recipe function f(sim::Simulation,ψ::Array,ind::Int,by::Union{Symbol,Function})
	markershape --> :rect
	markerstrokealpha --> 0
	seriestype --> :scatter
	q = quantile(abs.(ψ[:,ind]),PLOT_QUANTILE)
	if by ∈ [:real, :Real, :re, :Re, real]
		@series begin
			title --> "Re psi"
			markercolor --> :diverging
			colorbar --> false
			clim --> (-q,q).*PLOT_SCALE_FUDGE
			z = real(ψ[:,ind])
			marker_z --> z
			(sim.x[sim.interior],sim.y[sim.interior])
		end
	elseif by ∈ [:imag, :Imag, :im, :Im, imag]
		@series begin
			title --> "Im psi"
			markercolor --> :diverging
			colorbar --> false
			clim --> (-q,q).*PLOT_SCALE_FUDGE
			z = imag(ψ[:,ind])
			marker_z --> z
			(sim.x[sim.interior],sim.y[sim.interior])
		end
	elseif by ∈ [:abs2, :Abs2, abs2]
		@series begin
			title --> "|psi|²"
			markercolor --> :sequential
			colorbar --> false
			clim --> (0,q^2).*PLOT_SCALE_FUDGE
			z = abs2.(ψ[:,ind])
			marker_z --> z
			(sim.x[sim.interior],sim.y[sim.interior])
		end
	elseif by ∈ [:abs, :Abs, abs]
		@series begin
            title --> "|psi|"
			markercolor --> :sequential
			colorbar --> false
			z = abs.(ψ[:,ind])
			clim --> (0,q).*PLOT_SCALE_FUDGE
			marker_z --> z
			(sim.x[sim.interior],sim.y[sim.interior])
		end
	end
end

end # module

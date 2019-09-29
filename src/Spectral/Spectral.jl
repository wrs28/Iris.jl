"""
	module Spectral
"""
module Spectral

export eig_kl
export eig_knl
export eig_cf

using ..Common
using ArnoldiMethod
using ArnoldiMethodTransformations: AbstractSolver, PSOLVER, USOLVER, MSOLVER
# using Distributed
using LinearAlgebra
using NonlinearEigenproblems
# using ProgressMeter
# using Random
# using RecipesBase
using SparseArrays
# using Statistics

struct IrisLinSolverCreator <: LinSolverCreator end
struct IrisLinSolver{S,L,T} <: LinSolver
    LU::L
    x::Vector{T}
end
function NonlinearEigenproblems.create_linsolver(::IrisLinSolverCreator,nep,λ)
    M = compute_Mder(nep,λ)
    type = eltype(M)
    x = Vector{type}(undef,size(M,1))
    temp2 = Vector{type}(undef,size(M,1))
    S,LU = ArnoldiMethodTransformations.initialize_according_to_package(LUPACK[],issymmetric(M),type,M,x,temp2)
    return IrisLinSolver{S,typeof(LU),eltype(x)}(LU,x)
end
function NonlinearEigenproblems.lin_solve(solver::IrisLinSolver{S},b::Vector;tol=eps()) where S
    if S<:PSolver
        pardiso(solver.LU,solver.x,spzeros(eltype(b),length(b),length(b)),b)
    else
        ldiv!(solver.x, solver.LU, b)
    end
    return solver.x
end

include("1D/Spectral.jl")
# include("2D/Spectral.jl")
# include("3D/Spectral.jl")

# @recipe function f(sim::Simulation,ψ::Array,inds=1:size(ψ,2); by=nothing, structure=false)
# 	@assert issubset(inds,1:size(ψ,2)) "indices $inds inconsistent with size(ψ,2)=$(size(ψ,2))"
# 	legend --> false
# 	aspect_ratio --> 1
# 	n = length(inds)
#     if isnothing(by)
#         if structure
#             layout --> (1+n,3)
# 			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/sqrt(9+n^2+1)
#         else
#             layout --> (n,3)
# 			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/sqrt(9+n^2)
#         end
#     else
#         if structure
#             layout --> (3+n)
# 			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/sqrt(9+n^2)
#         else
#             layout --> (n)
# 			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/n
#         end
#     end
# 	if structure
# 		if isnothing(by)
# 			@series (sim,by)
# 		else
# 			@series (sim,real)
# 		end
# 	end
# 	for i ∈ 1:n
# 		if isnothing(by)
# 			bys = [:real, :imag, :abs2]
# 			for j ∈ eachindex(bys)
# 				if structure
# 					@series begin
# 						subplot --> 3i+j
# 						(sim,ψ,inds[i],bys[j])
# 					end
# 				else
# 					@series begin
# 						subplot --> 3(i-1)+j
# 						(sim,ψ,inds[i],bys[j])
# 					end
# 				end
# 			end
# 		else
# 			@series begin
# 				subplot --> i
# 				colorbar --> false
# 				(sim,ψ,inds[i],by)
# 			end
# 		end
# 	end
# end
# @recipe function f(sim::Simulation,ψ::Array,ind::Int,by::Union{Symbol,Function})
# 	markershape --> :rect
# 	markerstrokealpha --> 0
# 	seriestype --> :scatter
# 	q = quantile(abs.(ψ[:,ind]),PLOT_QUANTILE)
# 	if by ∈ [:real, :Real, :re, :Re, real]
# 		@series begin
# 			title --> "Re psi"
# 			markercolor --> :diverging
# 			colorbar --> false
# 			clim --> (-q,q).*PLOT_SCALE_FUDGE
# 			z = real(ψ[:,ind])
# 			marker_z --> z
# 			(sim.x[sim.interior],sim.y[sim.interior])
# 		end
# 	elseif by ∈ [:imag, :Imag, :im, :Im, imag]
# 		@series begin
# 			title --> "Im psi"
# 			markercolor --> :diverging
# 			colorbar --> false
# 			clim --> (-q,q).*PLOT_SCALE_FUDGE
# 			z = imag(ψ[:,ind])
# 			marker_z --> z
# 			(sim.x[sim.interior],sim.y[sim.interior])
# 		end
# 	elseif by ∈ [:abs2, :Abs2, abs2]
# 		@series begin
# 			title --> "|psi|²"
# 			markercolor --> :sequential
# 			colorbar --> false
# 			clim --> (0,q^2).*PLOT_SCALE_FUDGE
# 			z = abs2.(ψ[:,ind])
# 			marker_z --> z
# 			(sim.x[sim.interior],sim.y[sim.interior])
# 		end
# 	elseif by ∈ [:abs, :Abs, abs]
# 		@series begin
#             title --> "|psi|"
# 			markercolor --> :sequential
# 			colorbar --> false
# 			z = abs.(ψ[:,ind])
# 			clim --> (0,q).*PLOT_SCALE_FUDGE
# 			marker_z --> z
# 			(sim.x[sim.interior],sim.y[sim.interior])
# 		end
# 	end
# end

end # module

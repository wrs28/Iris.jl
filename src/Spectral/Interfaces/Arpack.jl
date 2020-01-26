module ArpackInterface

using Arpack
using LinearAlgebra
using SparseArrays

using ..Common
import ..AbstractEigenproblem
import ..AbstractLinearEigenproblem
import ..AbstractCFEigenproblem
import ..AbstractNonlinearEigenproblem
import ..HelmholtzProblem
import ..MaxwellProblem
import ..DEFAULT_LINEAR_EIGENSOLVER
import ..INDEX_OFFSET
import ..HelmholtzLEP
import ..HelmholtzCF
import ..MaxwellLEP
import ..MaxwellCF
import ..DOC_LEP_H
import ..DOC_LEP_M
import ..DOC_CF_H
import ..DOC_CF_M
import ..helmholtzeigen
import ..maxwelleigen
import ..orthogonalize!

################################################################################
# LEP

"""
    eigs(lep::AbstractLinearEigenproblem, ω, [ωs, ψs; nev=6, ncv=max(20,2*nev+1), which=:LM, tol=0.0, maxiter=300, sigma=ω², ritzvec=true, v0=zeros((0,))])

if frequencies `ωs::Vector` and fields `ψs::VectorField` are provided, susceptability is *saturated*
"""
function Arpack.eigs(lep::AbstractLinearEigenproblem, ω::Number, args...; kwargs...)
    A, B = lep(ω, args...)
    ω², ψ, nconv, niter, nmult, resid = eigs(spdiagm(0=>1 ./diag(B))*A; kwargs..., sigma=ω^2)
	if typeof(lep) <: HelmholtzLEP
		Ψ = ScalarField(lep,ψ)
	elseif typeof(lep) <: MaxwellLEP
		Ψ = ElectricField(lep,ψ)
	end
    return sqrt.(ω²), Ψ, nconv, niter, nmult, resid
end


################################################################################
# CF

"""
    eigs(cf::AbstractCFEigenproblem, ω, [ωs, ψs; η=0, nev=6, ncv=max(20,2*nev+1), which=:LM, tol=0.0, maxiter=300, sigma=η, ritzvec=true, v0=zeros((0,))])

if frequencies `ωs::Vector` and fields `ψs::ElectricField` are provided, susceptability is *saturated*
"""
function Arpack.eigs(cf::AbstractCFEigenproblem, ω::Number, args...; η::Number=0, kwargs...)
    A, B = cf(ω, args...)
    ηs, u, nconv, niter, nmult, resid = eigs(A, B; kwargs..., sigma=η)
	if typeof(cf) <: HelmholtzCF
		U = ScalarField(cf, u)
	elseif typeof(cf) <: MaxwellCF
		U = ElectricField(cf, u)
	end
    return ηs, U, nconv, niter, nmult, resid
end


################################################################################
# maxwelleigen and helmholtzeigen

function iris_eigen_arpack(
            lep::AbstractLinearEigenproblem,
            ω::Number,
			args...;
            verbose::Bool=false,
            kwargs...)

	ωs, ψs, nconv, niter, nmult, resid = eigs(lep, ω, args...; kwargs...)
    length(ωs) ≤ nconv || @warn "$(length(ωs) - nconv) evecs did not converge"
    normalize!(lep.simulation, ψs.values, lep.αεpFχ, size(ψs,1)÷2+INDEX_OFFSET) # Normalize according to (ψ₁,ψ₂)=δ₁₂
	if typeof(lep) <: HelmholtzLEP
		orthogonalize!(ψs,lep.simulation, ωs, lep.αεpFχ, 0, 0)
	elseif typeof(lep) <: MaxwellLEP
		orthogonalize!(ψs,lep.simulation, ωs, lep.αεpFχ, lep.ky, lep.kz)
		if all(iszero,(lep.ky,lep.kz)) normalize!(lep.simulation, ψs.values, lep.αεpFχ, size(ψs,1)÷2+INDEX_OFFSET) end
	end
	return ωs, ψs
end

if DEFAULT_LINEAR_EIGENSOLVER == :Arpack
	@doc """
	$DOC_LEP_H
	`nev` Number of eigenvalues (`6`);
	`v0` Starting vector (`zeros((0,))`);
	`maxiter` Maximum iterations (`300`);
	`ncv` Number of Krylov vectors (`max(20,2*nev+1)`);
	`tol` Tolerance is max of ε and `tol` (`0.0`);
	"""
	helmholtzeigen(lep::HelmholtzLEP, args...;kwargs...) = iris_eigen_arpack(lep,args...;kwargs...)

	@doc """
	$DOC_LEP_M
	`nev` Number of eigenvalues (`6`);
	`v0` Starting vector (`zeros((0,))`);
	`maxiter` Maximum iterations (`300`);
	`ncv` Number of Krylov vectors (`max(20,2*nev+1)`);
	`tol` Tolerance is max of ε and `tol` (`0.0`);
	"""
	maxwelleigen(lep::MaxwellLEP, args...;kwargs...) = iris_eigen_arpack(lep,args...;kwargs...)

end


function iris_eigen_arpack(
            cf::AbstractCFEigenproblem,
            ω::Number,
			args...;
			η::Number = 0,
            verbose::Bool=false,
            kwargs...)

	ηs, us, nconv, niter, nmult, resid = eigs(cf, ω, args...; kwargs...)
    length(ηs) ≤ nconv || @warn "$(length(ηs) - nconv) evecs did not converge"
    normalize!(cf.simulation, us, cf.F) # Normalize according to (ψ₁,ψ₂)=δ₁₂
	orthogonalize!(us,cf.simulation, ηs, cf.F, cf.ky, cf.kz)
	if all(iszero,(cf.ky,cf.kz)) normalize!(cf.simulation, us, cf.F) end
	return ηs, us
end

if DEFAULT_LINEAR_EIGENSOLVER == :Arpack
	@doc """
	$DOC_CF_H
	`nev` Number of eigenvalues (`6`);
	`v0` Starting vector (`zeros((0,))`);
	`maxiter` Maximum iterations (`300`);
	`ncv` Number of Krylov vectors (`max(20,2*nev+1)`);
	`tol` Tolerance is max of ε and `tol` (`0.0`);
	""" ->
	helmholtzeigen(cf::HelmholtzCF, args...;kwargs...) = iris_eigen_arpack(cf,args...;kwargs...)

	@doc """
	$DOC_CF_M
	`nev` Number of eigenvalues (`6`);
	`v0` Starting vector (`zeros((0,))`);
	`maxiter` Maximum iterations (`300`);
	`ncv` Number of Krylov vectors (`max(20,2*nev+1)`);
	`tol` Tolerance is max of ε and `tol` (`0.0`);
	""" ->
	maxwelleigen(cf::MaxwellCF, args...;kwargs...) = iris_eigen_arpack(cf,args...;kwargs...)

end

end #module

using .ArpackInterface

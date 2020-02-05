"""
Interface between `Iris` and `Arpack`
"""
module ArpackInterface

using Arpack
using ..Common

################################################################################
# LEP

using LinearAlgebra
using SparseArrays
import ..AbstractLinearEigenproblem
import ..HelmholtzLEP
import ..MaxwellLEP

"""
    eigs(lep, ω, [ωs, ψs; nev=6, ncv=max(20,2*nev+1), which=:LM, tol=0.0, maxiter=300, sigma=ω², ritzvec=true, v0=zeros((0,))]) -> (ω, [ψ,], nconv, niter, nmult, resid)

`lep` is an `AbstractLinearEigenproblem`, constructed from `HelmholtzLEP` or `MaxwellLEP`.
If frequencies `ωs::Vector` and fields `ψs::VectorField` are provided, susceptability is *saturated*.
Linear eigenproblem evaluates any boundary matching and dispersive domains at `ω`.
"""
function Arpack.eigs(
			lep::AbstractLinearEigenproblem,
			ω::Number,
			args...;
			kwargs...)

	# get A, B matrices from lep
    A, B = lep(ω, args...)
	# solve B\A for eigenpair
    ω²::Vector{ComplexF64}, ψ::Matrix{ComplexF64}, nconv, niter, nmult, resid = eigs(spdiagm(0=>1 ./diag(B))*A; kwargs..., sigma=ω^2)
	# wrap ψ as a VectorField, depending on Helmhotz vs Maxwell
	if typeof(lep) <: HelmholtzLEP
		Ψ = ScalarField(lep,ψ)
	elseif typeof(lep) <: MaxwellLEP
		Ψ = ElectricField(lep,ψ)
	end
	# return ω (eval is ω²)
    return sqrt.(ω²), Ψ, nconv, niter, nmult, resid
end


################################################################################
# CF

using LinearAlgebra
using SparseArrays
import ..AbstractCFEigenproblem
import ..HelmholtzCF
import ..MaxwellCF

"""
    eigs(cf, ω, [ωs, ψs; η=0, nev=6, ncv=max(20,2*nev+1), which=:LM, tol=0.0, maxiter=300, sigma=η, ritzvec=true, v0=zeros((0,))]) -> (η, [u,], nconv, niter, nmult, resid)

`cf` is an `AbstractCFEigenproblem`, constructed from `HelmholtzCF` or `MaxwellCF`.
If frequencies `ωs::Vector` and fields `ψs::ElectricField` are provided, susceptability is *saturated*.
CF eigenproblem evaluates any boundary matching and dispersive domains at `ω`.
"""
function Arpack.eigs(
			cf::AbstractCFEigenproblem,
			ω::Number,
			args...;
			η::Number=0,
			kwargs...)

	# get A, B from cf
    A, B = cf(ω, args...)
	# note Arpack can't handle B which is not positive semidefinite
		lmul!(-1,A)
		lmul!(-1,B)
		n = length(cf.simulation)
		sgn = sparse(I*1,n,n)
		rows = rowvals(B)
		vals = nonzeros(B)
		@inbounds for col ∈ 1:n
			@inbounds for j ∈ nzrange(B, col)
				row = rows[j]
				sgn[row,col] = flipsign(sgn[row,col],B.nzval[j])
				B.nzval[j] = abs(B.nzval[j])
			end
		end
	# now do actual solve
    ηs::Vector{ComplexF64}, u::Matrix{ComplexF64}, nconv, niter, nmult, resid = eigs(sgn*A, B; kwargs..., sigma=η)
	# wrap u as a VectorField, depending on Helmhotz vs Maxwell
	if typeof(cf) <: HelmholtzCF
		U = ScalarField(cf, u)
	elseif typeof(cf) <: MaxwellCF
		U = ElectricField(cf, u)
	end
    return ηs, U, nconv, niter, nmult, resid
end


################################################################################
# maxwelleigen and helmholtzeigen, LEP

using LinearAlgebra
import ..orthogonalize!
import ..INDEX_OFFSET

# applies when Arpack is default interface
function iris_eigen_arpack(
            lep::AbstractLinearEigenproblem,
            ω::Number,
			args...;
            verbose::Bool=false,
            kwargs...)

	# use previously defined eigs
	ωs, ψs, nconv, niter, nmult, resid = eigs(lep, ω, args...; kwargs...)
    length(ωs) ≤ nconv || @warn "$(length(ωs) - nconv) evecs did not converge"
	# Normalize according to (ψ₁,ψ₂)=δ₁₂
    normalize!(lep.simulation, ψs.values, lep.αεpFχ, size(ψs,1)÷2+INDEX_OFFSET)
	# Orthogonalize degenerate evecs
	if typeof(lep) <: HelmholtzLEP
		orthogonalize!(ψs,lep.simulation, ωs, lep.αεpFχ, 0, 0)
	elseif typeof(lep) <: MaxwellLEP
		orthogonalize!(ψs,lep.simulation, ωs, lep.αεpFχ, lep.ky, lep.kz)
		# for Maxwell ky=kz=0, there is inherent two-fold polarization degeneracy, this is cheap way to treat
		if all(iszero,(lep.ky,lep.kz)) normalize!(lep.simulation, ψs.values, lep.αεpFχ, size(ψs,1)÷2+INDEX_OFFSET) end
	end
	return ωs, ψs
end


import ..DEFAULT_LINEAR_EIGENSOLVER
import ..DOC_LEP_H
import ..DOC_LEP_M
import ..helmholtzeigen
import ..maxwelleigen

if DEFAULT_LINEAR_EIGENSOLVER == :Arpack
	@doc """
	$DOC_LEP_H
	`nev` Number of eigenvalues (`6`) |
	`v0` Starting vector (`zeros((0,))`) |
	`maxiter` Maximum iterations (`300`) |
	`ncv` Number of Krylov vectors (`max(20,2*nev+1)`) |
	`tol` Tolerance is max of ε and `tol` (`0.0`)
	"""
	helmholtzeigen(lep::HelmholtzLEP, args...;kwargs...) = iris_eigen_arpack(lep,args...;kwargs...)

	@doc """
	$DOC_LEP_M
	`nev` Number of eigenvalues (`6`) |
	`v0` Starting vector (`zeros((0,))`) |
	`maxiter` Maximum iterations (`300`) |
	`ncv` Number of Krylov vectors (`max(20,2*nev+1)`) |
	`tol` Tolerance is max of ε and `tol` (`0.0`)
	"""
	maxwelleigen(lep::MaxwellLEP, args...;kwargs...) = iris_eigen_arpack(lep,args...;kwargs...)
end

################################################################################
# maxwelleigen and helmholtzeigen, CF

function iris_eigen_arpack(
            cf::AbstractCFEigenproblem,
            ω::Number,
			args...;
			η::Number=0,
            verbose::Bool=false,
            kwargs...)

	# use previously defined eigs
	ηs, us, nconv, niter, nmult, resid = eigs(cf, ω, args...; η=η, kwargs...)
    length(ηs) ≤ nconv || @warn "$(length(ηs) - nconv) evecs did not converge"
	# Normalize according to (ψ₁,ψ₂)=δ₁₂ via CF.F, which may be saturated
    normalize!(cf.simulation, us, cf.F)
	# Orthogonalize degenerate evecs
	orthogonalize!(us,cf.simulation, ηs, cf.F, cf.ky, cf.kz)
	# for Maxwell ky=kz=0, there is inherent two-fold polarization degeneracy, this is cheap way to treat
	if all(iszero,(cf.ky,cf.kz)) normalize!(cf.simulation, us, cf.F) end
	return ηs, us
end

import ..DOC_CF_H
import ..DOC_CF_M

if DEFAULT_LINEAR_EIGENSOLVER == :Arpack
	@doc """
	$DOC_CF_H
	`nev` Number of eigenvalues (`6`) |
	`v0` Starting vector (`zeros((0,))`) |
	`maxiter` Maximum iterations (`300`) |
	`ncv` Number of Krylov vectors (`max(20,2*nev+1)`) |
	`tol` Tolerance is max of ε and `tol` (`0.0`)
	""" ->
	helmholtzeigen(cf::HelmholtzCF, args...;kwargs...) = iris_eigen_arpack(cf,args...;kwargs...)

	@doc """
	$DOC_CF_M
	`nev` Number of eigenvalues (`6`) |
	`v0` Starting vector (`zeros((0,))`) |
	`maxiter` Maximum iterations (`300`) |
	`ncv` Number of Krylov vectors (`max(20,2*nev+1)`) |
	`tol` Tolerance is max of ε and `tol` (`0.0`)
	""" ->
	maxwelleigen(cf::MaxwellCF, args...;kwargs...) = iris_eigen_arpack(cf,args...;kwargs...)
end

end #module

using .ArpackInterface

using ArnoldiMethod
using ArnoldiMethodTransformations

################################################################################
#LEP

"""
	partialschur(::AbstractLinearEigenproblem, ω; lupack=$DEFAULT_LUPACK, kwargs...) -> decomp, history

	partialschur(::AbstractLinearEigenproblem, ω, ωs::Array, ψs::ElectricField; lupack=$DEFAULT_LUPACK, kwargs...) -> decomp, history
"""
function ArnoldiMethod.partialschur(
			lep::AbstractLinearEigenproblem,
			ω::Number,
			args...;
			lupack::AbstractLUPACK=DEFAULT_LUPACK,
			kwargs...)

	A, B = lep(ω,args...)
	return partialschur(A, B, ω^2; diag_inv_B=true, lupack=lupack, kwargs...)
end

"""
	partialeigen(lep::AbstractLinearEigenproblem, decomp, ω) -> ωs, electricfields

order of `decomp` and `ω` does not matter
"""
function ArnoldiMethod.partialeigen(
			lep::AbstractLinearEigenproblem,
			ω::Number,
			decomp)

	ωs::Vector{ComplexF64}, ψ::Matrix{ComplexF64} = partialeigen(decomp, ω^2)
	normalize!(lep.sim, ψ, lep.εpχ, size(ψ,1)÷2+INDEX_OFFSET) # Normalize according to (ψ₁,ψ₂)=δ₁₂
	if typeof(lep) <: HelmholtzLEP
		Ψ = ScalarField(lep,ψ)
	elseif typeof(lep) <: MaxwellLEP
		Ψ = ElectricField(lep,ψ)
	end
	return sqrt.(ωs), Ψ
end
ArnoldiMethod.partialeigen(lep::AbstractLinearEigenproblem,decomp,ω::Number) = partialeigen(lep,ω,decomp)

function iris_eigen_arnoldimethod(
            lep::AbstractLinearEigenproblem,
            ω::Number,
			args...;
            verbose::Bool=false,
            lupack::AbstractLUPACK=DEFAULT_LUPACK,
            kwargs...)

	decomp, history = partialschur(lep, ω, args...; diag_inv_B=true, lupack=lupack, kwargs...)
    history.converged || @warn "$(history.nev - history.nconverged) evecs did not converge"
    verbose ? println(history) : nothing
    return partialeigen(lep, ω, decomp)
end

if DEFAULT_LINEAR_EIGENSOLVER == :ArnoldiMethod
@doc """
$doc_lep_h
`nev` Number of eigenvalues (`6`);
`verbose` Print convergence info (`false`)

Advanced keywords:
`lupack` LU solver (`$DEFAULT_LUPACK`);
`which` See `ArnoldiMethod` documentation (`LM()`);
`tol` Tolerance for convergence (`√eps`);
`mindim` Minimum Krylov dimension (`max(10, nev)`);
`maxdim` Maximum Krylov dimension (`max(20, 2nev)`);
`restarts` Maximum number of restarts (`200`)
""" ->
helmholtzeigen(lep::HelmholtzLEP,args...;kwargs...) = iris_eigen_arnoldimethod(lep,args...;kwargs...)

@doc """
$doc_lep_m
`nev` Number of eigenvalues (`6`);
`verbose` Print convergence info (`false`)

Advanced keywords:
`lupack` LU solver (`$DEFAULT_LUPACK`);
`which` See `ArnoldiMethod` documentation (`LM()`);
`tol` Tolerance for convergence (`√eps`);
`mindim` Minimum Krylov dimension (`max(10, nev)`);
`maxdim` Maximum Krylov dimension (`max(20, 2nev)`);
`restarts` Maximum number of restarts (`200`)
""" ->
maxwelleigen(lep::MaxwellLEP,args...;kwargs...) = iris_eigen_arnoldimethod(lep,args...;kwargs...)
end

################################################################################
#CF

"""
	partialschur(::AbstractCFEigenproblem, ω; η=0, lupack=$DEFAULT_LUPACK, kwargs...) -> decomp, history

	partialschur(::AbstractCFEigenproblem, ω, ωs, ψs; η=0, lupack=$DEFAULT_LUPACK, kwargs...) -> decomp, history
"""
function ArnoldiMethod.partialschur(
			cf::AbstractCFEigenproblem,
			ω::Number,
			args...;
			η::Number=0,
			lupack::AbstractLUPACK=DEFAULT_LUPACK,
			kwargs...)

	A, B = cf(ω,args...)
	return partialschur(A, B, η; diag_inv_B=true, lupack=lupack, kwargs...)
end

"""
	partialeigen(cf::AbstractCFEigenproblem, decomp, η) -> ηs, us
"""
function ArnoldiMethod.partialeigen(
			cf::AbstractCFEigenproblem,
			η::Number,
			decomp)

	ηs::Vector{ComplexF64}, us::Matrix{ComplexF64} = partialeigen(decomp, η)
	normalize!(cf.maxwell.sim, us, cf.F)
	if typeof(cf) <: HelmholtzCF
		U = ScalarField(cf,us)
	elseif typeof(cf) <: MaxwellCF
		U = ElectricField(cf,us)
	end
	return sqrt.(ηs), U
end
ArnoldiMethod.partialeigen(cf::AbstractCFEigenproblem,decomp,η::Number) = partialeigen(cf,η,decomp)

function iris_eigen_arnoldimethod(
            cf::AbstractCFEigenproblem,
            ω::Number,
			args...;
            η::Number=0,
            verbose::Bool=false,
            lupack::AbstractLUPACK=DEFAULT_LUPACK,
			kwargs...)

    decomp, history = partialschur(cf, ω, args...; η=η, diag_inv_B=false, lupack=lupack, kwargs...)
    history.converged || @warn "$(history.nev-history.nconverged) eigenvectors failed to converge"
    verbose ? println(history) : nothing
    return partialeigen(cf, η, decomp)
end

if DEFAULT_LINEAR_EIGENSOLVER == :ArnoldiMethod
@doc """
$doc_cf_h
`nev` Number of eigenvalues (`6`);
`verbose` Print convergence info (`false`);

Advanced keywords:
`lupack` LU solver (`$DEFAULT_LUPACK`);
`which` See `ArnoldiMethod` documentation (`LM()`);
`tol` Tolerance for convergence (`√eps`);
`mindim` Minimum Krylov dimension (`max(10, nev)`);
`maxdim` Maximum Krylov dimension (`max(20, 2nev)`);
`restarts` Maximum number of restarts (`200`)
""" ->
helmholtzeigen(cf::HelmholtzCF,args...;kwargs...) = iris_eigen_arnoldimethod(cf,args...;kwargs...)

@doc """
$doc_cf_m
`nev` Number of eigenvalues (`6`);
`verbose` Print convergence info (`false`);

Advanced keywords:
`lupack` LU solver (`$DEFAULT_LUPACK`);
`which` See `ArnoldiMethod` documentation (`LM()`);
`tol` Tolerance for convergence (`√eps`);
`mindim` Minimum Krylov dimension (`max(10, nev)`);
`maxdim` Maximum Krylov dimension (`max(20, 2nev)`);
`restarts` Maximum number of restarts (`200`)
""" ->
maxwelleigen(cf::MaxwellCF,args...;kwargs...) = iris_eigen_arnoldimethod(cf,args...;kwargs...)
end

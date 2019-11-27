using ArnoldiMethod
using ArnoldiMethodTransformations

################################################################################
#LEP

"""
	partialschur(::MaxwellLEP, ω; lupack=$DEFAULT_LUPACK, kwargs...) -> decomp, history

	partialschur(::MaxwellLEP, ω, ωs::Array, ψs::ElectricField; lupack=$DEFAULT_LUPACK, kwargs...) -> decomp, history
"""
function ArnoldiMethod.partialschur(
			lep::MaxwellLEP,
			ω::Number,
			args...;
			lupack::AbstractLUPACK=DEFAULT_LUPACK,
			kwargs...)

	A, B = lep(ω,args...)
	return partialschur(A, B, ω^2; diag_inv_B=true, lupack=lupack, kwargs...)
end

"""
	partialeigen(lep::MaxwellLEP, decomp, ω) -> ωs, electricfields

order of `decomp` and `ω` does not matter
"""
function ArnoldiMethod.partialeigen(
			lep::MaxwellLEP,
			ω::Number,
			decomp)

	ωs::Vector{ComplexF64}, ψ::Matrix{ComplexF64} = partialeigen(decomp, ω^2)
	normalize!(lep.sim, ψ, lep.εpχ, size(ψ,1)÷2+INDEX_OFFSET) # Normalize according to (ψ₁,ψ₂)=δ₁₂
	return sqrt.(ωs), ElectricField(lep,ψ)
end
ArnoldiMethod.partialeigen(lep::MaxwellLEP,decomp,ω::Number) = partialeigen(lep,ω,decomp)

function maxwelleigen_arnoldimethod(
            lep::MaxwellLEP,
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
$(doc_lep)
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
maxwelleigen(lep::MaxwellLEP,args...;kwargs...) = maxwelleigen_arnoldimethod(lep,args...;kwargs...)
end

################################################################################
#CF

"""
	partialschur(::MaxwellCF, ω; η=0, lupack=$DEFAULT_LUPACK, kwargs...) -> decomp, history

	partialschur(::MaxwellCF, ω, ωs, ψs; η=0, lupack=$DEFAULT_LUPACK, kwargs...) -> decomp, history
"""
function ArnoldiMethod.partialschur(
			cf::MaxwellCF,
			ω::Number,
			args...;
			η::Number=0,
			lupack::AbstractLUPACK=DEFAULT_LUPACK,
			kwargs...)

	A, B = cf(ω,args...)
	return partialschur(A, B, η; diag_inv_B=true, lupack=lupack, kwargs...)
end

"""
	partialeigen(cf::MaxwellCF, decomp, ω; η=0) -> ωs, electricfields
"""
function ArnoldiMethod.partialeigen(
			cf::MaxwellCF,
			ω::Number,
			decomp;
			η::Number=0)

	ηs::Vector{ComplexF64}, us::Matrix{ComplexF64} = partialeigen(decomp, η)
	normalize!(cf.maxwell.sim, us, cf.F)
	return sqrt.(ηs), ElectricField(cf,us)
end
ArnoldiMethod.partialeigen(cf::MaxwellCF,decomp,ω::Number;kwargs...) = partialeigen(cf,ω,decomp;kwargs...)

function maxwelleigen_arnoldimethod(
            cf::MaxwellCF,
            ω::Number,
			args...;
            η::Number=0,
            verbose::Bool=false,
            lupack::AbstractLUPACK=DEFAULT_LUPACK,
			kwargs...)

    decomp, history = partialschur(cf, ω, args...; η=η, diag_inv_B=false, lupack=lupack, kwargs...)
    history.converged || @warn "$(history.nev-history.nconverged) eigenvectors failed to converge"
    verbose ? println(history) : nothing
    return partialeigen(cf, ω, decomp; η=η)
end

if DEFAULT_LINEAR_EIGENSOLVER == :ArnoldiMethod
@doc """
$(doc_cf)
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
maxwelleigen(cf::MaxwellCF,args...;kwargs...) = maxwelleigen_arnoldimethod(cf,args...;kwargs...)
end

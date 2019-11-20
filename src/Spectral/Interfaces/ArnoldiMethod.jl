using ArnoldiMethod

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
"""
function ArnoldiMethod.partialeigen(
			lep::MaxwellLEP,
			ω::Number,
			decomp)

	ωs::Vector{ComplexF64}, ψ::Matrix{ComplexF64} = partialeigen(decomp, ω^2)
	normalize!(lep.sim, ψ, lep.εpχ, size(ψ,1)÷2+INDEX_OFFSET) # Normalize according to (ψ₁,ψ₂)=δ₁₂
	return sqrt.(ωs), ElectricField(lep,ψ)
end

# for now maxwell_eigen(::MaxwellLEP, ...) uses ArnoldiMethod
"""
	maxwelleigen_arnoldi(lep::MaxwellLEP,ω,[ωs,ψs; verbose=false, lupack=$DEFAULT_LUPACK, kwargs...)
"""
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


################################################################################
#CF


"""
	partialschur(::MaxwellCF, ω; lupack=$DEFAULT_LUPACK, kwargs...) -> decomp, history
	partialschur(::MaxwellCF, ω, ωs::Array, ψs::ElectricField; lupack=$DEFAULT_LUPACK, kwargs...) -> decomp, history
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
	partialeigen(cf::MaxwellCF, decomp, ω) -> ωs, electricfields
"""
function ArnoldiMethod.partialeigen(
			cf::MaxwellCF,
			ω::Number,
			decomp)

	ηs::Vector{ComplexF64}, us::Matrix{ComplexF64} = partialeigen(decomp, ω^2)
	normalize!(cf.maxwell.sim, us, cf.F)
	return sqrt.(ηs), ElectricField(cf,us)
end

"""
	maxwelleigen_arnoldimethod(cf::MaxwellCF,ω,[ωs,ψs; η=0, verbose=false, lupack=$DEFAULT_LUPACK, kwargs...)
"""
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
    return partialeigen(decomp,η)
end

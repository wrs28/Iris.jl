using Arpack

"""
    eigs(lep::MaxwellLEP, ω, [ωs, ψs; nev=6, ncv=max(20,2*nev+1), which=:LM, tol=0.0, maxiter=300, sigma=ω², ritzvec=true, v0=zeros((0,))])

if frequency vector `ωs` and electric fields `ψs` are provided, F is *saturated* F
"""
function Arpack.eigs(lep::MaxwellLEP, ω::Number, args...; kwargs...)
    A, B = lep(ω, args...)
    ω², ψ, nconv, niter, nmult, resid = eigs(spdiagm(0=>1 ./diag(B))*A; kwargs..., sigma=ω^2)
    return sqrt.(ω²), ElectricField(lep.simulation,ψ), nconv, niter, nmult, resid
end

"""
    eigs(cf::MaxwellCF, ω, [ωs, ψs; η=0, nev=6, ncv=max(20,2*nev+1), which=:LM, tol=0.0, maxiter=300, sigma=η, ritzvec=true, v0=zeros((0,))])

if frequency vector `ωs` and electric fields `ψs` are provided, F is *saturated* F
"""
function Arpack.eigs(cf::MaxwellCF, ω::Number, args...; η::Number=0, kwargs...)
    A, B = cf(ω, args...)
    ηs, u, nconv, niter, nmult, resid = eigs(A, B; kwargs..., sigma=η)
    return ηs, ElectricField(cf.simulation, u), nconv, niter, nmult, resid
end


"""
	maxwelleigen_arpack(lep::MaxwellLEP,ω,[ωs,ψs; verbose=false, kwargs...)
"""
function maxwelleigen_arpack(
            lep::MaxwellLEP,
            ω::Number,
			args...;
            verbose::Bool=false,
            kwargs...)

	ωs, ψs, nconv, niter, nmult, resid = eigs(lep, ω, args...; kwargs...)
    length(ωs) ≤ nconv || @warn "$(length(ωs) - nconv) evecs did not converge"
    normalize!(lep.sim, ψs.values, lep.αεpχ, size(ψs,1)÷2+INDEX_OFFSET) # Normalize according to (ψ₁,ψ₂)=δ₁₂
	return ωs, ψs
end


"""
	maxwelleigen_arpack(lep::MaxwellCF,ω,[ωs,ψs; verbose=false, kwargs...)
"""
function maxwelleigen_arpack(
            cf::MaxwellCF,
            ω::Number,
			args...;
			η::Number = 0,
            verbose::Bool=false,
            kwargs...)

	ηs, us, nconv, niter, nmult, resid = eigs(cf, ω, args...; kwargs...)
    length(ηs) ≤ nconv || @warn "$(length(ηs) - nconv) evecs did not converge"
    normalize!(cf.sim, us.values, cf.F, size(us,1)÷2+INDEX_OFFSET) # Normalize according to (ψ₁,ψ₂)=δ₁₂
	return ηs, us
end

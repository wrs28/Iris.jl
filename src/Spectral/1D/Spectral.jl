function LinearAlgebra.normalize!(sim::Simulation{1},ψ,B)
    n = length(sim.domain_indices)
    lattices = map(d->d.lattice,sim.domains)
    dx = map(i->lattices[i].dx,sim.domain_indices)
    for i ∈ 1:size(ψ,2)
        𝒩² = zero(eltype(ψ))
        for j ∈ 1:size(ψ,1)
            𝒩² += (ψ[j,i].^2)*B[j,j]*dx[mod1(j,n)]
        end
        ψ[:,i] /= sqrt(𝒩²)*exp(complex(0,angle(ψ[end÷2-1,i])))
    end
    return nothing
end


"""
    eig_kl(sim, ω; ky=0, kz=0, verbose=false, kwargs...]) -> ω, ψ

Linear eigenvalue solve for `ω`.

Keyword `verbose` defaults to `false`.
See docs for ArnoldiMethodWrapper for details of `kwargs`.
"""
function eig_kl(sim::Simulation{1,Symmetric,T},
            ω::Number;
            ky::Number=0,
            kz::Number=0,
            verbose::Bool=false,
            kwargs...
            ) where T

    A,B = maxwell_lep(sim;ky=ky,kz=kz)(ω)
    decomp, history = partialschur(A, B, ω^2; diag_inv_B=true, kwargs...)
    history.converged || @warn "$(history.nev - history.nconverged) did not converge"
    verbose ? println(history) : nothing

    λ::Array{T,1}, v::Array{T,2} = partialeigen(decomp, ω^2)
    normalize!(sim,v,B) # Normalize according to (ψ₁,ψ₂)=δ₁₂
    return sqrt.(λ), v
end


"""
    eig_knl(sim, k, ka=0, kb=0; method=contour_beyn, nk=3, display=false, quad_n=100, kwargs...) -> k,ψ
"""
function eig_knl(
            sim::Simulation{1,Symmetric,T},
            ω::Number;
            ky::Number=0,
            kz::Number=0,
            verbose::Bool = false,
            check_consistency::Bool = true,
            Schur_fact::Bool = false,
            align_sparsity_patterns::Bool = false,
            nlmethod::Function = iar,
            logger::Integer = Int(verbose),
            kwargs...
            ) where T

    spmf_kwargs = (
        :check_consistency=>check_consistency,
        :Schur_fact=>Schur_fact,
        :align_sparsity_patterns=>align_sparsity_patterns)
    nep = maxwell_nep(sim;
        ky=ky,
        kz=kz,
        :check_consistency => check_consistency,
        :Schur_fact => Schur_fact,
        :align_sparsity_patterns => align_sparsity_patterns)

    λ::Array{T,1}, v::Array{T,2} = nlmethod(nep; σ=ω, kwargs...,logger=logger)
    return λ, v
end


"""
    eig_cf(sim, k, [ka=0, kb=0; η=0, verbose, lupack, kwargs...]) -> η, u

Linear CF eigenvalues closest to `η`

Keyword `verbose` defaults to `false`.
Keyword `lupack` defaults to `:auto` and contrls which package is used in the eigensolver.
See docs for ArnoldiMethodWrapper for details.
"""
function eig_cf(sim::Simulation{1,Symmetric,T},
            ω::Number;
            ky::Number=0,
            kz::Number=0,
            η::Number=0,
            verbose::Bool=false,
            kwargs...
            ) where T

    A,B = maxwell_lep(sim)(ω)
    A = A-B*ω^2
    F = -spdiagm(0=>vcat(sim.α[1].*sim.F,sim.α[1].*sim.F,sim.α[1].*sim.F))
    decomp, history = partialschur(A, F*ω^2, η; diag_inv_B=false, kwargs...)
    history.converged || @warn "$(history.nev-history.nconverged) eigenvectors failed to converge"
    verbose ? println(history) : nothing

    λ::Array{T,1}, v::Array{T,2} = partialeigen(decomp,η)
    normalize!(sim,v,B)
    return λ, v
end

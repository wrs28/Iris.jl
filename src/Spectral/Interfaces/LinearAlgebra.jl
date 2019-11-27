using LinearAlgebra

function LinearAlgebra.normalize!(sim::Simulation{1},ψ,B)
    @inbounds for j ∈ 1:size(ψ,2)
        𝒩² = zero(eltype(ψ))
        # @fastmath @inbounds @simd
        for i ∈ 1:size(ψ,1)
            𝒩² += (ψ[i,j]^2)*B[i,i]*sim.dx
        end
        ψ[:,j] /= sqrt(𝒩²)
    end
    return nothing
end

function LinearAlgebra.normalize!(sim::Simulation{1},ψ,B,ind::Int)
    normalize!(sim,ψ,B)
    @fastmath @inbounds @simd for i ∈ 1:size(ψ,2) ψ[:,i] /= cis(angle(ψ[ind,i])) end
    return nothing
end

LinearAlgebra.eigen(prob::AbstractMaxwellEigenproblem, args...; kwargs...) =
    maxwelleigen(prob, args...; kwargs...)

using LinearAlgebra


function LinearAlgebra.normalize!(sim::Simulation{1},ψ,B)
    n = length(sim.domain_indices)
    lattices = map(d->d.lattice,sim.domains)
    dx = map(i->lattices[i].dx,sim.domain_indices)
    for i ∈ 1:size(ψ,2)
        𝒩² = zero(eltype(ψ))
        for j ∈ 1:size(ψ,1)
            𝒩² += (ψ[j,i].^2)*B[j,j]*dx[mod1(j,n)]
        end
        ψ[:,i] /= sqrt(𝒩²)
    end
    return nothing
end

function LinearAlgebra.normalize!(sim::Simulation{1},ψ,B,ind::Int)
    n = length(sim.domain_indices)
    normalize!(sim,ψ,B)
    for i ∈ 1:size(ψ,2) ψ[:,i] /= cis(angle(ψ[ind,i])) end
    return nothing
end


LinearAlgebra.eigen(prob::T, args...; kwargs...)  where T<:Union{MaxwellLEP, MaxwellCF, MaxwellNEP} =
    maxwelleigen(prob, args...; kwargs...)

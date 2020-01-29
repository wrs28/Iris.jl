module LinearAlgebraInterface

export _measure

using LinearAlgebra
using ..Common

import ..AbstractEigenproblem
import ..AbstractLinearEigenproblem
import ..AbstractCFEigenproblem
import ..AbstractNonlinearEigenproblem
import ..DEFAULT_LINEAR_EIGENSOLVER
import ..HelmholtzLEP
import ..MaxwellLEP

IrisSymmetric = Common.Symmetric

"""
    normalize!(simulation, ψ, B, [index])

normalize electric field `ψ` with respect to inner product `B`. If `index` is provided,[]
fixes `ψ[index,:]` to be real.
"""
function LinearAlgebra.normalize!(sim::Simulation,ψ,B)
    @inbounds for j ∈ 1:size(ψ,2)
        𝒩² = zero(eltype(ψ))
        # @fastmath @inbounds @simd
        for i ∈ 1:size(ψ,1)
            𝒩² += (ψ[i,j]^2)*B[i,i]*_measure(sim)
        end
        ψ[:,j] /= sqrt(𝒩²)
    end
    return nothing
end

function LinearAlgebra.normalize!(sim::Simulation,ψ,B,ind::Int)
    normalize!(sim,ψ,B)
    @fastmath @inbounds @simd for i ∈ 1:size(ψ,2) ψ[:,i] /= cis(angle(ψ[ind,i])) end
    return nothing
end

"""
    eigen(::AbstractEigenproblem, args...; kwargs...)

See [`helmholtzeigen`](@ref), [`maxwelleigen`](@ref)
"""
function LinearAlgebra.eigen(prob::AbstractEigenproblem, args...; kwargs...)
    if problemtype(prob) <: HelmholtzProblem
        return helmholtzeigen(prob, args...; kwargs...)
    elseif problemtype(prob) <: MaxwellProblem
        return maxwelleigen(prob, args...; kwargs...)
    end
end

function _measure(sim::Simulation{1})
    if typeof(sim.lattice)<:Lattice{1,Cartesian}
        D = sim.dx
    end
    return D
end

function _measure(sim::Simulation{2,IrisSymmetric})
    if typeof(sim.lattice)<:Lattice{2,Cartesian}
        D = sim.dx*sim.dy
    end
    return D
end

end # module

using .LinearAlgebraInterface

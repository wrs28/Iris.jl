module DivGradCurls

export Gradient
export Divergence
export Curl

using ..Common
# import ..Common.Symmetric
using LinearAlgebra
using SparseArrays

struct Gradient{N}
    components::NTuple{N,SparseMatrixCSC{Float64,Int}}
end

Gradient(components::Vararg{T,N}) where {T<:AbstractArray,N} = Gradient{N}(components)

function Gradient(sim::Simulation{1,Common.Symmetric})
    dx = sim.dx
    n = length(sim)
    m = n+1
    ∇ = spdiagm(m,n,0=>vcat(2/dx,fill(1/dx,n-1)),-1=>vcat(fill(-1/dx,n-1),-2/dx))
    return Gradient(∇)
end

function Gradient(sim::Simulation{2,Common.Symmetric})
    dx = sim.dx
    dy = sim.dy
    n = length(sim)
    m = n+1
    ∇ = spdiagm(m,n,0=>vcat(2/dx,fill(1/dx,n-1)),-1=>vcat(fill(-1/dx,n-1),-2/dx))
    return Gradient(∇)
end

function LinearAlgebra.mul!(Y::VectorField{N,N}, A::Gradient{N}, B::ScalarField{N}) where N
    for i ∈ 1:N
        mul!(component(i,Y), A.components[i], B.values)
    end
    return Y
end

struct Divergence{N}
    gradient::Gradient{N}
end

function Divergence(sim::Simulation{1,Common.Symmetric})
    dx = sim.dx
    m = length(sim)
    n = length(sim)+1
    ∇ = spdiagm(m,n,0=>vcat(fill(-1/dx,m)),1=>vcat(fill(1/dx,m)))
    return Divergence(Gradient(∇))
end

function LinearAlgebra.mul!(Y::ScalarField{N}, A::Divergence{N}, B::VectorField{N,N}) where N
    for i ∈ 1:N
        mul!(Y.values, A.gradient.components[i], component(i,B))
    end
    return Y
end

struct Curl{N}

end

function Curl(sim::Simulation{N}) where N

end

end # module

using .DivGradCurls

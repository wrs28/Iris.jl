module DivGradCurls

export Gradient
export Divergence
export Curl

using ..Common
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
    nx, ny = sim.latticesize
    mx = nx+1; my = ny+1
    ∇ₓ = spdiagm(mx,nx,0=>vcat(2/dx,fill(1/dx,nx-1)),-1=>vcat(fill(-1/dx,nx-1),-2/dx))
    ∇ᵤ = spdiagm(my,ny,0=>vcat(2/dy,fill(1/dy,ny-1)),-1=>vcat(fill(-1/dy,ny-1),-2/dy))
    return Gradient(kron(sparse(I,ny,ny),∇ₓ),kron(∇ᵤ,sparse(I,nx,nx)))
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

function Divergence(sim::Simulation{2,Common.Symmetric})
    dx = sim.dx
    dy = sim.dy
    mx, my = sim.latticesize
    nx = mx+1
    ny = my+1
    ∇ₓ = spdiagm(mx,nx,0=>vcat(fill(-1/dx,mx)),1=>vcat(fill(1/dx,mx)))
    ∇ᵤ = spdiagm(my,ny,0=>vcat(fill(-1/dy,my)),1=>vcat(fill(1/dy,my)))
    return Divergence(Gradient(kron(sparse(I,my,my),∇ₓ),kron(∇ᵤ,sparse(I,mx,mx))))
end


struct Curl{N}

end

function Curl(sim::Simulation{N}) where N

end

################################################################################
# extend mul!

function LinearAlgebra.mul!(Y::NTuple{N,ScalarField{N}}, A::Gradient{N}, B::ScalarField{N}) where N
    for i ∈ 1:N
        mul!(Y[i], A.components[i], B.values)
    end
    return Y
end

function LinearAlgebra.mul!(Y::ScalarField{N}, A::Divergence{N}, B::NTuple{N,ScalarField{N}}) where N
    if N==1
        mul!(Y.values, A.gradient.components[1], B[1].values)
    else
        fill!(Y.values,0)
        for i ∈ 1:N
            mul!(Y.values, A.gradient.components[i], B[i],1,1)
        end
    end
    return Y
end

end # module

using .DivGradCurls

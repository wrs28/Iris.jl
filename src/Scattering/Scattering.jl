module Scattering

export solve_scatter,
green_fn


using ..IrosBase
using ..Spectral
using ArnoldiMethodTransformations
using LinearAlgebra
using SparseArrays


function solve_scatter(sim::Simulation,k::Number,j::AbstractArray)
    L = Laplacian(sim,k,0,0)
    A = L.D + L.M*k^2
    return solve_scatter_core(A,j)
end


function green_fn(sim::Simulation,k::Number,x::Real,y::Real)
    inds = IrosBase.Tessellations.get_surrounding_sites(sim.tessellation,x,y)
    weights = IrosBase.Simulations.generate_weights(x,y,sim.x[inds],sim.y[inds])
    j = zeros(sim.n)
    x = sim.x[inds]; y = sim.y[inds]
    v1 = [x[2]-x[1],y[2]-y[1]]
    v2 = [x[3]-x[2],y[3]-y[2]]
    cosθ = (norm(v1+v2)^2 - norm(v1)^2 - norm(v2)^2)/2/norm(v1)/norm(v2)
    sinθ = sqrt(1-cosθ^2)
    area = norm(v1)*norm(v2)*sinθ/2
    j[inds] = weights./area
    return solve_scatter(sim,k,j)
end


function solve_scatter(sim::Simulation,kx::Number,ky::Number,shape::AbstractShape,f::Function)
    k = hypot(kx,ky)
    x = sim.x[sim.interior]
    y = sim.y[sim.interior]
    L = Laplacian(sim,k,0,0)
    φ = (f.(kx,ky,x,y).*shape.(x,y))
    j = (L.D+sparse(I*k^2,sim.n,sim.n))*φ
    A = L.D + L.M*k^2
    scattered_field = solve_scatter_core(A,-j)
    total_field = scattered_field + φ
    return total_field, scattered_field, φ
end
solve_scatter(sim::Simulation,k::Vector,shape::AbstractShape,f::Function) = solve_scatter(sim,k[1],k[2],shape,f)


function solve_scatter_core(A::SparseMatrixCSC,j::AbstractArray)
    A_lu! = ArnoldiMethodWrapper.ShiftAndInvert(A,0)
    y = Array{ComplexF64}(undef,length(j))
    solve_scatter_core!(y,A,j)
    return y
end
function solve_scatter_core!(y::AbstractArray,A::SparseMatrixCSC,j::AbstractArray)
    @assert length(j)==size(A,1) "length(j)=$(length(j)) must equal size of system $(size(A,1))"
    A_lu! = ArnoldiMethodWrapper.ShiftAndInvert(A,0)
    A_lu!(y,j)
    return nothing
end

end # module

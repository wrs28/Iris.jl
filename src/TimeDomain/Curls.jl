module Curls

export Curl

using ...Common
using SparseArrays

struct Curl{N,T}
	curlB::SparseMatrixCSC{T,Int}
	curlE::SparseMatrixCSC{T,Int}
end

function Curl(sim::Simulation{1})

	∂ₓE = E_derivative(length(sim),sim.dx)
	posE = map(p->p.x,sim.x)
	PE = invperm(sortperm(posE))

	∂ₓB = B_derivative(length(sim),sim.dx)
	posB = expand_B_pos(posE,sim.dx)
	PB = invperm(sortperm(posB))

	permute!(∂ₓE,PB,PE) # not sure why this is the correct permutation
	permute!(∂ₓB,PE,PB) # but it works

	I3 = sparse([3,2],[2,3],[1,-1],3,3)
	curlB = kron(I3,∂ₓB)
	curlE = kron(I3,∂ₓE)

	return Curl{1,Float64}(curlB,curlE)
end

function Curl(sim::Simulation{2}) end

function Curl(sim::Simulation{3}) end

function B_derivative(n::Int,dx::Float64)
	rowsm = 1:n
	rowsp = 1:n
	colsm = 1:n
	colsp = 2:n+1
	vals = ones(n)
	return sparse(vcat(rowsm,rowsp),vcat(colsm,colsp),vcat(-vals,vals)/dx,n,n+1)
end

function E_derivative(n::Int,dx::Float64)
	rowsm = 2:n
	rowsp = 2:n
	colsm = 1:n-1
	colsp = 2:n
	vals = ones(n-1)
	return sparse(vcat([1,n+1],rowsm,rowsp),vcat([1,n],colsm,colsp),vcat([2,-2],-vals,vals)/dx,n+1,n)
end

expand_B_pos(sim::Simulation{1}) = expand_B_pos(map(p->p.x,sim.x),sim.dx)
expand_B_pos(posE::Vector,dx::Float64) = vcat(posE.-dx/2,maximum(posE)+dx/2)

end # module

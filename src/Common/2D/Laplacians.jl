"""
    module Laplacians
"""
module Laplacians

export Laplacian

using ...Defaults
import ..Shapes
using ..Simulations
using SparseArrays


"""
	struct Laplacian

	Laplacian(sim,k) --> laplacian

fields:
	`D`: Laplacian operator with PMLs (total kinetic energy)
	`K`: part of the Laplacian due to ordinary coupling between sites
	`Σ`: part of the Laplacian due to the self-energy from the boundary conditions
	`M`: mass matrix (scaled dielectric)
	`Mf`: auxilliary mass matrix (scaled pump)
	`scaling`: scalings associated with PMLs and/or polar coordinates
	`k`: the frequency which was used to construct the PML's.
	`ka`, `kb`: the Bloch vectors used to construct the Laplacian.
"""
struct Laplacian{TK,TKX,TKY}
	D::SparseMatrixCSC{ComplexF64,Int64} # total kinetic operator
	K::SparseMatrixCSC{ComplexF64,Int64} # coupling term
	Σ::SparseMatrixCSC{ComplexF64,Int64} # self-energy
	M::SparseMatrixCSC{ComplexF64,Int64} # mass matrix
	Mf::SparseMatrixCSC{ComplexF64,Int64} # mass matrix
	scaling::SparseMatrixCSC{ComplexF64,Int64}
	k::Float64
	ka::TKX
	kb::TKY

	function Laplacian(sim::Simulation,k::Number,kx::Number=0,ky::Number=0)
		k = real(k)
		sx,sy = boundary_layer_site(sim)
		HX = spdiagm(0=>map(z->(1+1im*z/k),sx))
		HY = spdiagm(0=>map(z->(1+1im*z/k),sy))

		Kx,Σx = laplacian_x(sim,k)
		Ky,Σy = laplacian_y(sim,k)
		K = HY*Kx+HX*Ky
		Σ = HY*Σx+HX*Σy
		D = K + Σ
		scaling = HX*HY
		M = scaling*spdiagm(0=>sim.ε[sim.interior])
		Mf = scaling*spdiagm(0=>sim.F[sim.interior])
		return new{typeof(k),typeof(kx),typeof(ky)}(D, K, Σ, M, Mf, scaling, k, kx, ky)
	end

	function Base.show(io::IO,lap::Laplacian)
		println(io,"Laplacian:")
		println(io,"\toptimal frequency: k=$(lap.k)")
		println(io,"\tBloch phases: ka=$(lap.ka), kb=$(lap.ka)")
		println(io,"\tD (kinetic): $(nnz(lap.D)) nonzeros")
		println(io,"\tK (coupling): $(nnz(lap.K)) nonzeros")
		println(io,"\tΣ (self-energy): $(nnz(lap.Σ)) nonzeros")
		println(io,"\tM, Mf (mass): $(nnz(lap.M)) nonzeros")
	end
end
# r = lds.r0 + ij[i][1]*lds.dr
# r_half₋ = lds.r0 + (ij[i][1]-1/2)*lds.dr
# r_half₊ = lds.r0 + (ij[i][1]+1/2)*lds.dr

function laplacian_x(sim::Simulation,k::Real)
	N = sum(sim.interior)

	link_weight = complex.(sim.link_weight[sim.link_x_bool])
	s = boundary_layer_link_x(sim)
	h = @. 1/(1+1im*s/k)
	for i ∈ eachindex(link_weight)
		if !iszero(s[i])
			link_weight[i] *= h[i]
		end
	end
	start = sim.link_start[sim.link_x_bool]
	stop = sim.link_stop[sim.link_x_bool]
	rows = vcat(start, start)
	cols = vcat( stop, start)
	vals = vcat(link_weight,-link_weight)
	Kx = sparse(rows,cols,vals,N,N,+)


	self_weight = complex.(sim.self_weight[sim.self_x_bool])
	s = boundary_layer_self_x(sim)
	h = @. 1/(1+1im*s/k)
	for i ∈ eachindex(self_weight)
		if !iszero(s[i])
			self_weight[i] *= h[i]
		end
	end
	self = sim.self_index[sim.self_x_bool]
	Σx = sparse(self,self,self_weight,N,N,+)
	return Kx,Σx
end


function laplacian_y(sim::Simulation,k::Real)
	N = sum(sim.interior)

	link_weight = complex.(sim.link_weight[sim.link_y_bool])
	s = boundary_layer_link_y(sim)
	h = @. 1/(1+1im*s/k)
	for i ∈ eachindex(link_weight)
		if !iszero(s[i])
			link_weight[i] *= h[i]
		end
	end
	start = sim.link_start[sim.link_y_bool]
	stop = sim.link_stop[sim.link_y_bool]
	rows = vcat(start, start)
	cols = vcat( stop, start)
	vals = vcat(link_weight,-link_weight)
	Ky = sparse(rows,cols,vals,N,N,+)

	self_weight = complex.(sim.self_weight[sim.self_y_bool])
	s = boundary_layer_self_y(sim)
	h = @. 1/(1+1im*s/k)
	for i ∈ eachindex(self_weight)
		if !iszero(s[i])
			self_weight[i] *= h[i]
		end
	end
	self = sim.self_index[sim.self_y_bool]
	Σy = sparse(self,self,self_weight,N,N,+)
	return Ky,Σy
end


function boundary_layer_site(sim::Simulation)
	domains = sim.domains
	x = sim.x[sim.interior]
	y = sim.y[sim.interior]
	domain_index = sim.domain_index[sim.interior]
	return boundary_layer(domains,x,y,domain_index)
end

boundary_layer_link_x(sim::Simulation) = boundary_layer_link(sim,sim.link_x_bool)[1]
boundary_layer_link_y(sim::Simulation) = boundary_layer_link(sim,sim.link_y_bool)[2]
function boundary_layer_link(sim::Simulation,link_half_bool)
	domains = sim.domains
	x = sim.link_half_x[link_half_bool]
	y = sim.link_half_y[link_half_bool]
	domain_index = sim.domain_index[sim.link_start[link_half_bool]]
	return boundary_layer(domains,x,y,domain_index)
end


boundary_layer_self_x(sim::Simulation) = boundary_layer_self(sim,sim.self_x_bool)[1]
boundary_layer_self_y(sim::Simulation) = boundary_layer_self(sim,sim.self_y_bool)[2]
function boundary_layer_self(sim::Simulation,self_half_bool)
	domains = sim.domains
	x = sim.self_half_x[self_half_bool]
	y = sim.self_half_y[self_half_bool]
	domain_index = sim.domain_index[sim.self_index[self_half_bool]]
	return boundary_layer(domains,x,y,domain_index)
end

function boundary_layer(domains,x,y,domain_index)
	Σx = Array{ComplexF64}(undef,length(x))
	Σy = Array{ComplexF64}(undef,length(y))
	boundary_layer!(Σx,Σy,domains,x,y,domain_index)
	return Σx, Σy
end
function boundary_layer!(Σx,Σy,domains,x,y,domain_index)
	for i ∈ eachindex(x)
		Σx[i],Σy[i] = boundary_layer(domains[domain_index[i]],x[i],y[i])
	end
	return nothing
end
function boundary_layer(domain,x,y)
	x,y = Shapes.unrotate(x,y,domain.boundary.shape)
	N = Shapes.get_number_of_sides(domain.boundary.shape)
	if N==4
	    bl1 = domain.boundary.bls[1]
		bl2 = domain.boundary.bls[2]
		bl3 = domain.boundary.bls[3]
		bl4 = domain.boundary.bls[4]
		σx = bl1(x,y)+bl2(x,y)
		σy = bl3(x,y)+bl4(x,y)
	elseif N==2
		bl1 = domain.boundary.bls[1]
		bl2 = domain.boundary.bls[2]
		σx = bl1(x,y)+bl2(x,y)
        σy = complex(0.)
	elseif N==1
		bl1 = domain.boundary.bls[1]
		σx = bl1(x,y)
		σy = complex(0.)
	else
		σx, σy =complex(0.), complex(0.)
	end
    return σx,σy
end


end #module

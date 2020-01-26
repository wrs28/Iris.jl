module Laplacians2DSymmetric

using ..Lattices
using ..Points
using SparseArrays

import ..Symmetric, ..Unsymmetric
import LinearAlgebra: I
import ..Laplacian

################################################################################
# CARTESIAN
"""
	Laplacian(::Lattice{2}, α::Vector, α_half::Vector) -> lap
"""
function Laplacian{Symmetric}(
			lattice::Lattice{2,Cartesian},
			αx::Vector{ComplexF64},
			αy::Vector{ComplexF64},
			α_halfx::Vector{ComplexF64},
			α_halfy::Vector{ComplexF64})

	dx = lattice.dx
	dy = lattice.dy

	Nx = length(αx)
	Ny = length(αy)
	N = Nx*Ny

	α_halfx⁻¹ = 1 ./α_halfx
	α_halfy⁻¹ = 1 ./α_halfy

	# note that x_half[1] = x[1] - dx/2, so that D₋½=D(x_half[1])
	# and that D½=D(x_half[2])
	rows0 = 1:Nx	; cols0 = 1:Nx	; vals0 = - α_halfx⁻¹[rows0] - α_halfx⁻¹[rows0.+1]
	rowsm = 2:Nx	; colsm = 1:Nx-1; valsm = α_halfx⁻¹[1 .+ colsm]
	rowsp = 1:Nx-1	; colsp = 2:Nx	; valsp = α_halfx⁻¹[colsp]
	∂ₓα⁻¹∂ₓ = sparse(vcat(rowsm,rows0,rowsp),vcat(colsm,cols0,colsp),vcat(valsm,vals0,valsp)/lattice.dx^2,Nx,Nx)

	rows0 = 1:Ny	; cols0 = 1:Ny	; vals0 = - α_halfy⁻¹[rows0] - α_halfy⁻¹[rows0.+1]
	rowsm = 2:Ny	; colsm = 1:Ny-1; valsm = α_halfy⁻¹[1 .+ colsm]
	rowsp = 1:Ny-1	; colsp = 2:Ny	; valsp = α_halfy⁻¹[colsp]
	∂ᵤα⁻¹∂ᵤ = sparse(vcat(rowsm,rows0,rowsp),vcat(colsm,cols0,colsp),vcat(valsm,vals0,valsp)/lattice.dy^2,Ny,Ny)

	Ax = spdiagm(0=>αx)
	Ay = spdiagm(0=>αy)

	∇² = kron(Ay, ∂ₓα⁻¹∂ₓ) + kron(∂ᵤα⁻¹∂ᵤ, Ax)

	return Laplacian{2,Symmetric}(∇², kron(Ay,Ax))
end

################################################################################
# POLAR

function Laplacian{Symmetric}(lattice::Lattice{2,Polar}, α::Vector{ComplexF64}, α_half::Vector{ComplexF64})#, nnm, nnp, indices, interior, surface)
	throw("not implemented for 2D Polar yet")
	N = length(α)

	dr = lattice.dr
	r = (1/2 .+ (0:N-1))*dr

	α_half⁻¹ = 1 ./α_half
	r_half = (0:N)*dr

	∂ₓα⁻¹∂ₓ_bulk = _bulk_laplacian(lattice,N,interior,surface,nnm,nnp,α_half⁻¹)
	∂ₓα⁻¹∂ₓ_surface = _surface_laplacian(lattice,N,interior,surface,nnm,nnp,α_half⁻¹,indices)

	∂ₓ_bulk = _bulk_derivative(lattice,N,interior,surface,nnm,nnp)
	∂ₓ_surface = _surface_derivative(lattice,N,interior,surface,nnm,nnp,indices)
	∂ₓ = ∂ₓ_bulk + ∂ₓ_surface
	i∂ₓ = 1im*∂ₓ

	return Laplacian{2}(∂ₓα⁻¹∂ₓ)
end

function _bulk_laplacian(lattice::Lattice{2},N,interior,surface,nnm,nnp,α_half⁻¹)
	rows = findall(interior .& .!surface)
	colsm = nnm[1][interior .& .!surface]
	colsp = nnp[1][interior .& .!surface]
	valsm = α_half⁻¹[rows] # note that x_half[1] = x[1] - dx/2, so that D₋½=D(x_half[1])
	valsp = α_half⁻¹[colsp] #  and D½=D(x_half[2])
	cols0 = rows
	vals0 = - valsm - valsp
	return sparse(vcat(rows,rows,rows),vcat(colsm,cols0,colsp),vcat(valsm,vals0,valsp)/lattice.dx^2,N,N)
end
function _surface_laplacian(lattice::Lattice{2},N,interior,surface,nnm,nnp,α_half⁻¹,indices)
	rows_sf = findall(surface .& interior)
	cols_sfm = map(i->findfirst(isequal(indices[i]-CartesianIndex(1,)),indices),rows_sf)
	rows_sfm = [rows_sf[findfirst(.!isnothing.(cols_sfm))]]
	cols_sfm = [cols_sfm[findfirst(.!isnothing.(cols_sfm))]]
	vals_sfm = α_half⁻¹[rows_sfm] # note that x_half[1] = x[1] - dx/2, so that D₋½=D(x_half[1])
	cols_sfp = map(i->findfirst(isequal(indices[i]+CartesianIndex(1,)),indices),rows_sf)
	rows_sfp = [rows_sf[findfirst(.!isnothing.(cols_sfp))]]
	cols_sfp = [cols_sfp[findfirst(.!isnothing.(cols_sfp))]]
	vals_sfp = α_half⁻¹[cols_sfp] # and D½=D(x_half[2])
	cols_sf0 = rows_sf
	vals_sf0 = - α_half⁻¹[[0,1] .+ rows_sf] - α_half⁻¹[[0,1] .+ rows_sf]
	return sparse(vcat(rows_sfm,rows_sf,rows_sfp),vcat(cols_sfm,cols_sf0,cols_sfp),vcat(vals_sfm,vals_sf0,vals_sfp)/lattice.dx^2,N,N)
end

function _bulk_derivative(lattice::Lattice{2},N,interior,surface,nnm,nnp)
	rows = findall(interior .& .!surface)
	colsm = nnm[1][interior .& .!surface]
	valsm = -ones(size(colsm))
	colsp = nnp[1][interior .& .!surface]
	valsp = ones(size(colsm))
	return sparse(vcat(rows,rows),vcat(colsm,colsp),vcat(valsm,valsp)/lattice.dx,N,N)
end
function _surface_derivative(lattice::Lattice{2},N,interior,surface,nnm,nnp,indices)
	rows_sf = findall(surface .& interior)
	cols_sfm = map(i->findfirst(isequal(indices[i]-CartesianIndex(1,)),indices),rows_sf)
	rows_sfm = [rows_sf[findfirst(.!isnothing.(cols_sfm))]]
	cols_sfm = [cols_sfm[findfirst(.!isnothing.(cols_sfm))]]
	vals_sfm = -ones(size(rows_sfm))
	cols_sfp = map(i->findfirst(isequal(indices[i]+CartesianIndex(1,)),indices),rows_sf)
	rows_sfp = [rows_sf[findfirst(.!isnothing.(cols_sfp))]]
	cols_sfp = [cols_sfp[findfirst(.!isnothing.(cols_sfp))]]
	vals_sfp = ones(size(rows_sfm))
	return sparse(vcat(rows_sfm,rows_sfp),vcat(cols_sfm,cols_sfp),vcat(vals_sfm,vals_sfp)/lattice.dx,N,N)
end


################################################################################
# SPHERICAL
Laplacian{Symmetric}(lattice::Lattice{2,Spherical}, α::Vector{ComplexF64}, α_half::Vector{ComplexF64}, nnm, nnp, indices, interior, surface) = throw("not implemented for 2D Spherical yet")

end

using .Laplacians2DSymmetric

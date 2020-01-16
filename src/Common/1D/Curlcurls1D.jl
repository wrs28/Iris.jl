"""
	Curlcurl(::Lattice{1}, α::Vector, α_half::Vector, nnm, nnp, indices, interior, surface)
"""
function Curlcurl{Symmetric}(lattice::Lattice{1,Cartesian}, α::Vector{ComplexF64}, α_half::Vector{ComplexF64}, nnm, nnp, indices, interior, surface)

	dx = lattice.dx
	α_half⁻¹ = 1 ./α_half
	∂ₓα⁻¹∂ₓ = Laplacian{Symmetric}(lattice,α,α_half).l0

	N = length(α)
	∂ₓ_bulk = _bulk_derivative(lattice,N,interior,surface,nnm,nnp)
	∂ₓ_surface = _surface_derivative(lattice,N,interior,surface,nnm,nnp,indices)
	∂ₓ = ∂ₓ_bulk + ∂ₓ_surface
	i∂ₓ = 1im*∂ₓ

	return Curlcurl{1,Symmetric}(∂ₓα⁻¹∂ₓ,i∂ₓ,spdiagm(0=>α),spdiagm(0=>α))
end

"""
	(cc::Curlcurl{1})(ky, kz) = ∇×∇×
"""
function (cc::Curlcurl{1})(ky::Real,kz::Real)

	ky² = ky^2
	kz² = kz^2
	kykz = ky*kz

	cc0 = kron(sparse([2,3],[2,3],[-1,-1],3,3),cc.cc0)
	cc1 = kron(sparse([2,3,1,1],[1,1,2,3],[ky,kz,ky,kz],3,3),cc.cc1)
	cc2 = kron(sparse([1,2,3,2,3],[1,2,2,3,3],[ky²+kz²,kz²,-kykz,-kykz,ky²],3,3),cc.cc2)

	return cc0+cc1+cc2
end


function _bulk_laplacian(lattice::Lattice{1},N,interior,surface,nnm,nnp,α_half⁻¹)
	rows = findall(interior .& .!surface)
	colsm = nnm[1][interior .& .!surface]
	colsp = nnp[1][interior .& .!surface]
	valsm = α_half⁻¹[rows] # note that x_half[1] = x[1] - dx/2, so that D₋½=D(x_half[1])
	valsp = α_half⁻¹[colsp] #  and D½=D(x_half[2])
	cols0 = rows
	vals0 = - valsm - valsp
	return sparse(vcat(rows,rows,rows),vcat(colsm,cols0,colsp),vcat(valsm,vals0,valsp)/lattice.dx^2,N,N)
end
function _surface_laplacian(lattice,N,interior,surface,nnm,nnp,α_half⁻¹,indices)
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

function _bulk_derivative(lattice::Lattice{1},N,interior,surface,nnm,nnp)
	rows = findall(interior .& .!surface)
	colsm = nnm[1][interior .& .!surface]
	valsm = -ones(size(colsm))
	colsp = nnp[1][interior .& .!surface]
	valsp = ones(size(colsm))
	return sparse(vcat(rows,rows),vcat(colsm,colsp),vcat(valsm,valsp)/lattice.dx,N,N)
end
function _surface_derivative(lattice::Lattice{1},N,interior,surface,nnm,nnp,indices)
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







# function _domain_wall_laplacian(lattice::Lattice{1},N,interior,domain_wall,nnm,nnp,α_half⁻¹,indices)
# 	rows_dw = findall(domain_wall .& interior)
# 	cols_dwm = map(i->findfirst(isequal(indices[i]-CartesianIndex(1,)),indices),rows_dw)
# 	vals_dwm = α_half⁻¹[cols_dwm]
# 	cols_dwp = map(i->findfirst(isequal(indices[i]+CartesianIndex(1,)),indices),rows_dw)
# 	vals_dwp = α_half⁻¹[cols_dwp]
# 	cols_dw0 = rows_dw
# 	vals_dw0 = - α_half⁻¹[cols_dwm] - α_half⁻¹[cols_dwp]
# 	return sparse(vcat(rows_dw,rows_dw,rows_dw),vcat(cols_dwm,cols_dw0,cols_dwp),vcat(vals_dwm,vals_dw0,vals_dwp)/lattice.dx^2,N,N)
# end
# function _domain_wall_derivative(lattice::Lattice{1},N,interior,domain_wall,nnm,nnp,α_half⁻¹,indices)
# 	rows_dw = findall(domain_wall .& interior)
# 	cols_dwm = map(i->findfirst(isequal(indices[i]-CartesianIndex(1,)),indices),rows_dw)
# 	vals_dwm = -ones(size(cols_dwm))
# 	cols_dwp = map(i->findfirst(isequal(indices[i]+CartesianIndex(1,)),indices),rows_dw)
# 	vals_dwp = ones(size(cols_dwp))
# 	return sparse(vcat(rows_dw,rows_dw),vcat(cols_dwm,cols_dwp),vcat(vals_dwm,vals_dwp)/lattice.dx,N,N)
# end

import ..Curlcurls._bulk_derivative
import ..Curlcurls._domain_wall_derivative
import ..Curlcurls._surface_derivative

function Curl(lattice::Lattice{1},nnm,nnp,indices,interior,surface,domain_wall) where T
	dx = lattice.dx
	N = length(interior)

	∂ₓ_bulk = _bulk_derivative(lattice,N,interior,domain_wall,surface,nnm,nnp)
	∂ₓ_domain_wall = _domain_wall_derivative(lattice,N,interior,domain_wall,nnm,nnp,1,indices)
	∂ₓ_surface = _surface_derivative(lattice,N,interior,surface,nnm,nnp,indices)
	return Curl{1}(∂ₓ_bulk + ∂ₓ_domain_wall + ∂ₓ_surface)
end

(cc::Curl{1})() = kron(sparse([3,2],[2,3],[1,-1],3,3),cc.curl)

function (cc::Curl{1})(ky::Real,kz::Real)
	cc0 = cc()
	cc1 = kron(sparse([2,3,1,1],[1,1,2,3],1im*[kz,-ky,-kz,ky],3,3),sparse(I,size(cc0)...))
	return cc0+cc1
end

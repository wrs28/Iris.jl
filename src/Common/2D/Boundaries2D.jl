function _bls_by_shape(bls::AbstractBL{SIDE},c::Circle) where SIDE
	if SIDE == 1
		stop = c.radius
		start = stop + bls.depth
	else
		throw("circle only has side=1, evaluated at $side")
	end
	return bls(start=start,stop=stop)
end

# Square is just special rectangle
_bls_by_shape(bls::AbstractBL,s::Square) = _bls_by_shape(bls,Rectangle(s))

function _bls_by_shape(bls::AbstractBL{SIDE},shape::Rectangle) where SIDE
	if SIDE==1
		stop = shape.origin.x - shape.a/2
		start = stop + bls.depth
	elseif SIDE==2
		stop = shape.origin.x + shape.a/2
		start = stop - bls.depth
	elseif SIDE==3
		stop = shape.origin.y - shape.b/2
		start = stop + bls.depth
	elseif SIDE==4
		stop = shape.origin.y + shape.b/2
		start = stop - bls.depth
	end
	return bls(start=start, stop=stop)
end


################################################################################
# Plotting

@recipe function f(bnd::Boundary{2})
	bnd.shape
end

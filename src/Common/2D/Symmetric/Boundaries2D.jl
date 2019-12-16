function _bls_by_shape(bls::AbstractBL{SIDE},c::Circle) where SIDE
	if SIDE == 1
		stop = c.radius
		start = stop + bls.depth
	else
		throw("circle only has side=1, evaluated at $side")
	end
	return bls(start=start,stop=stop)
end

# function _bls_by_shape(bls::AbstractBL{SIDE},s::Square) where SIDE
# 	if SIDE == 1
# 		stop = s.a
# 		start = stop + bls.depth
# 	else
# 		throw("circle only has side=1, evaluated at $side")
# 	end
# 	return bls(start=start,stop=stop)
# end

_bls_by_shape(bls,s::Square) = _bls_by_shape(bls,Rectangle(s))

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

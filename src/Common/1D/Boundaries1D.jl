function _bls_by_shape(bls::AbstractBL{SIDE},i::Interval) where SIDE
	if SIDE == 1
		stop = i.start
		start = stop + bls.depth
	else
		stop = i.stop
		start = stop - bls.depth
	end
	return bls(start=start,stop=stop)
end

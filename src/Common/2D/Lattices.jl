function lattice_primitives(φ,constants::NTuple{2})
    sinφ, cosφ = sincos(φ)
    R = @SMatrix [  cosφ    -sinφ;
                    sinφ    cosφ]
    return (Point(R[:,1]),Point(R[:,2])), R
end

Base.getindex(lat::Lattice{2},ind1,ind2) = map(i,j->lat[i,j],ind1,ind2)
function Base.getindex(lat::Lattice{2},ind1::Real,ind2::Real)
    throw("not written yet")
end

latticeindex(lat::Lattice{2},p::Point{2}) = throw("not written yet")

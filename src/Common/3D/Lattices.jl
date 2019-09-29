function lattice_primitives(φ,θ,γ,constants::NTuple{3})
    sinφ, cosφ = sincos(φ)
    sinθ, cosθ = sincos(θ)
    sinγ, cosγ = sincos(γ)
    Rz = @SMatrix [ cosφ    -sinφ   0;
                    sinφ    cosφ    0;
                    0       0       1]
    Rx = @SMatrix [ 1       0       0;
                    0       cosθ    -sinθ;
                    0       sinθ    cosθ]
    Rz′ = @SMatrix [cosγ    -sinγ   0;
                    sinγ    cosγ    0;
                    0       0       1]
    R = Rz′*Rx*Rz
    return (Point(R[:,1]),Point(R[:,2]),Point(R[:,3])), R
end

Base.getindex(lat::Lattice{3},ind1,ind2,ind3) = map(i,j,k->lat[i,j,k],ind1,ind2,ind3)
function Base.getindex(lat::Lattice{2},ind1::Real,ind2::Real,ind3::Real)
    throw("not written yet")
end

latticeindex(lat::Lattice{3},p::Point{3}) = throw("not written yet")

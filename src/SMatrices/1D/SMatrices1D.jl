#TODO: Fix SMatrix for nonnormal incidence
"""
    smatrix(::Simulation{1}, ω; [ky=0, kz=0, lupack=$DEFAULT_LUPACK]) -> S
"""
function smatrix(sim::Simulation{1},
            ω::AbstractVector;
            ky::Real=0,
            kz::Real=0,
            lupack::AbstractLUPACK=DEFAULT_LUPACK)

    spc = SMatrixContainer(sim,ω,lupack,float(ky),float(kz))
    S = SharedArray{ComplexF64,3}((ndims(spc),ndims(spc),length(ω)); init=spc)
    return ScatteringMatrix{ndims(spc),typeof(ω)}(sdata(S),ω)
end

function SMatrixContainer(sim::Simulation{1},ω,lupack,ky,kz)
    N = sum(which_channels(sim))
    N>0 || throw("no scattering for closed system")
    return SMatrixContainer{N,typeof(sim),typeof(ω),typeof(lupack)}(sim,ω,lupack,ky,kz)
end

function which_channels(sim::Simulation{1})
    bc1,bc2 = sim.boundary.bcs
    bl1,bl2 = sim.boundary.bls
    bch1, bch2 = BCHermiticity(bc1), BCHermiticity(bc2)
    a = falses(2)
    if typeof(bl1)<:AbstractComplexBL || bch1<:NonHermitianBC
        a[1] = true
    end
    if typeof(bl2)<:AbstractComplexBL || bch2<:NonHermitianBC
        a[2] = true
    end
    return a
end

# One-sided, one-dimensional s-matrix core
@inline function dimensional_smatrix!(S::SharedArray, spc::SMatrixContainer{1,TSIM}, indices) where TSIM<:Simulation{1}
    ls = MaxwellLS(spc.sim, mean(spc.ω[indices]), 1, 0)
    a = which_channels(spc.sim)

    # @fastmath @inbounds @simd
    for i ∈ indices
        if a[1]
            update!(ls,1,0,spc.ω[i])
            alu = lu(ls.maxwell.A,spc.lupack)
            scattering!(ls,alu)
            S[1,1,i] = ls.solution.total.left[2]*sqrt(spc.ω[i])
        elseif a[2]
            update!(ls,0,1,spc.ω[i])
            alu = lu(ls.maxwell.A,spc.lupack)
            scattering!(ls,alu)
            S[1,1,i] = ls.solution.total.right[2]*sqrt(spc.ω[i])
        end
    end
    return nothing
end

# Two-sided, one-dimensional smatrix core
@inline function dimensional_smatrix!(S::SharedArray, spc::SMatrixContainer{2,TSIM}, indices) where TSIM<:Simulation{1}
    ls = MaxwellLS(spc.sim, mean(spc.ω[indices]), 1, 0)
    a = which_channels(spc.sim)

    # @fastmath @inbounds @simd
    for i ∈ indices
        update!(ls,1,0,spc.ω[i])
        alu = lu(ls.maxwell.A,spc.lupack)
        scattering!(ls,alu)
        S[1,1,i] = ls.solution.scattered.left[2]*sqrt(ls.equivalent_source.channelflux[1])
        S[1,2,i] = ls.solution.total.right[2]*sqrt(ls.equivalent_source.channelflux[2])

        update!(ls,0,1)
        scattering!(ls,alu)
        S[2,1,i] = ls.solution.total.left[2]*sqrt(ls.equivalent_source.channelflux[1])
        S[2,2,i] = ls.solution.scattered.right[2]*sqrt(ls.equivalent_source.channelflux[2])
    end
    return nothing
end

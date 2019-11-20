module SMatrices

export smatrix

using ..Common
using ..Scattering
using Distributed
using Interpolations
using LaTeXStrings
using LinearAlgebra
using RecipesBase
using SharedArrays

import ..Common.PRINTED_COLOR_LIGHT
import ..Common.PRINTED_COLOR_NUMBER
import ..Scattering.EquivalentSource
import ..Scattering.scattering!

struct ScatteringMatrix{N,TW}
    S::Array{ComplexF64,3}
    ω::TW
end

struct SMatrixContainer{TSIM,TW,TLU} <: Function
    sim::TSIM
    ω::TW
    lupack::TLU
    ky::Float64
    kz::Float64
end

"""
    smatrix(::Simulation, ω; [kx,ky,kz,lupack=$(DEFAULT_LUPACK)]) -> S

Frequency `ω` is a (real) number, vector, or range. Any argument order is valid.
Scattering matrix `S` is `N`-by-`N`-by-`M`, where `M=length(ω)` and `N` is the number
of channels.

Noteworthy: `S` is easily plotted with Plots.jl: `plot(S)`. Also automatically parallelizes over frequency `ω` if `nprocs()>1`.
"""
function smatrix(sim::Simulation{1},
            ω::AbstractVector;
            ky::Real=0,
            kz::Real=0,
            lupack::AbstractLUPACK=DEFAULT_LUPACK)

    spc = SMatrixContainer(sim,ω,lupack,float(ky),float(kz))
    S = SharedArray{ComplexF64,3}((2,2,length(ω)); init=spc)
    return ScatteringMatrix{2,typeof(ω)}(sdata(S),ω)
end
smatrix(ω, sim::Simulation; kwargs...) = smatrix(sim, ω; kwargs...)
smatrix(sim::Simulation, ω::Real; kwargs...) = smatrix(sim, [ω]; kwargs...)


@inline function (spc::SMatrixContainer)(S::SharedArray)
    if indexpids(S)==0
        return nothing
    else
        idx = indexpids(S)
        nchunks = length(procs(S))
        splits = [round(Int, s) for s in range(0, stop=size(S,3), length=nchunks+1)]
        indices = splits[idx]+1:splits[idx+1]
    end

    aL = [1,0]; aR = [0,1]; ω = spc.ω
    sol = ScatteringSolution(spc.sim, ω[1], aL)
    scL = EquivalentSource(spc.sim, sol.ω, aL)
    scR = EquivalentSource(spc.sim, sol.ω, aR)
    M = maxwell(spc.sim);

    @fastmath @inbounds @simd for i ∈ indices
        M(ω[i])
        alu = lu(M.A,spc.lupack)

        scL.ω = ω[i]; scR.ω = ω[i]
        jL = scL(M.D²,spc.ky,spc.kz)
        jR = scR(M.D²,spc.ky,spc.kz)

        scattering!(sol,M,scL,jL,spc.lupack,alu)
        S[1,1,i] = sol.sct.left[2]*sqrt(ω[i])
        S[1,2,i] = sol.tot.right[2]*sqrt(ω[i])

        scattering!(sol,M,scR,jR,spc.lupack,alu)
        S[2,1,i] = sol.tot.left[2]*sqrt(ω[i])
        S[2,2,i] = sol.sct.right[2]*sqrt(ω[i])
    end
    return nothing
end


# Pretty Printing
################################################################################
function Base.show(io::IO,s::ScatteringMatrix{N}) where N
    printstyled(io,"ScatteringMatrix ",color=PRINTED_COLOR_LIGHT)
    print(io,"($N-by-$N) @ $(length(s.ω)) frequencies")
end


# Plotting
################################################################################
# first smooth version when ω is LinRange
@recipe f(s::ScatteringMatrix{2,TW};n=1001) where TW<:LinRange = s,n
@recipe f(n::Integer,s::ScatteringMatrix{2,TW}) where TW<:LinRange = s,n
@recipe function f(s::ScatteringMatrix{2,TW},n::Integer) where TW<:LinRange
    ωs = LinRange(s.ω[1],s.ω[end],n)
    @series begin
        label := L"|t|^2"
        T = CubicSplineInterpolation(s.ω,abs2.(s.S[1,2,:]))
        ωs, T(ωs)
    end
    @series begin
        label := L"|r_L|^2"
        RL = CubicSplineInterpolation(s.ω,abs2.(s.S[1,1,:]))
        ωs, RL(ωs)
    end
    @series begin
        label := L"|r_R|^2"
        RR = CubicSplineInterpolation(s.ω,abs2.(s.S[2,2,:]))
        ωs, RR(ωs)
    end
end

# next jagged version when ω is Array
@recipe function f(s::ScatteringMatrix{2})
    @series begin
        label := L"|t|^2"
        s.ω, abs2.(s.S[1,2,:])
    end
    @series begin
        label := L"|r_L|^2"
        s.ω, abs2.(s.S[1,1,:])
    end
    @series begin
        label := L"|r_R|^2"
        s.ω, abs2.(s.S[2,2,:])
    end
end

end # module

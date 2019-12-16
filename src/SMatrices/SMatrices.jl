"""
Computation of Scattering Matrices

Accelerated if `nprocs()>1`.
"""
module SMatrices

export smatrix

dimension_files = (
    "1D/SMatrices1D.jl",
    )

using ..Common
import ..Common.AbstractComplexBL

using ..Scattering
using Distributed
using Interpolations
using LaTeXStrings
using LinearAlgebra
using RecipesBase
using SharedArrays
using Statistics

import ..Scattering.ScatteringSolution
import ..Scattering.EquivalentSource

"""
    ScatteringMatrix{N}

`N`×`N` Scattering Matrix container with fields `S` (scattering matrix), `ω` (frequencies)
"""
struct ScatteringMatrix{N,TW}
    S::Array{ComplexF64,3}
    ω::TW
end

"""
    SMatrixContainer{N}

Auxilliary for computation of `N`×`N` scattering matrix.
"""
struct SMatrixContainer{N,TSIM,TW,TLU} <: Function
    sim::TSIM
    ω::TW
    lupack::TLU
    ky::Float64
    kz::Float64
end

foreach(include,dimension_files)

Base.ndims(smc::SMatrixContainer{N}) where N = N

"""
    smatrix(::Simulation, ω; [kx, ky, kz, lupack=$(DEFAULT_LUPACK)]) -> S

Frequency `ω` is a (real) number, vector, or range.
Scattering matrix `S` is N-by-N-by-M, where M=`length(ω)` and N is the number
of channels.

Automatically parallelizes over frequency `ω` if `nprocs()>1`.

Noteworthy: `S` is easily plotted with Plots.jl: `plot(S)`.
"""
smatrix(sim::Simulation, ω::Real; kwargs...) = smatrix(sim, [ω]; kwargs...)

# core smatrix computation taking place during SharedMatrix initialization
function (spc::SMatrixContainer)(S::SharedArray)
    if indexpids(S)==0
        return nothing
    else
        idx = indexpids(S)
        nchunks = length(procs(S))
        splits = [round(Int, s) for s in range(0, stop=size(S,3), length=nchunks+1)]
        indices = splits[idx]+1:splits[idx+1]
    end
    dimensional_smatrix!(S,spc,indices)
    return nothing
end


################################################################################
# Pretty Printing
import ..Common.PRINTED_COLOR_LIGHT
import ..Common.PRINTED_COLOR_NUMBER

function Base.show(io::IO,s::ScatteringMatrix{N}) where N
    printstyled(io,"ScatteringMatrix ",color=PRINTED_COLOR_LIGHT)
    print(io,"($N-by-$N) @ $(length(s.ω)) frequencies")
end


################################################################################
# Plotting

# first smooth version when ω is LinRange
@recipe f(s::ScatteringMatrix{N,TW};n=1001) where {N,TW<:LinRange} = s,n
@recipe function f(s::ScatteringMatrix{1,TW},n::Integer) where TW<:LinRange
    ωs = LinRange(s.ω[1],s.ω[end],n)
    @series begin
        label := L"|r|^2"
        RL = CubicSplineInterpolation(s.ω,abs2.(s.S[1,1,:]))
        ωs, RL(ωs)
    end
end
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
@recipe function f(s::ScatteringMatrix{1,Vector})
    @series begin
        label := L"|r|^2"
        s.ω, abs2.(s.S[1,1,:])
    end
end
@recipe function f(s::ScatteringMatrix{2,Vector})
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

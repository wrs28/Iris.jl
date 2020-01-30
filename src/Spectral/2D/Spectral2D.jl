"""
Two-dimensional NEP
"""
module Spectral2D

using ..Common
using SparseArrays
using NonlinearEigenproblems

import ..HelmholtzNEP
import ..MaxwellNEP

"""
    HelmholtzNEP(sim) -> nep

`sim` is a `Simulation`
"""
function HelmholtzNEP(
            sim::Simulation{1,C,T};
            kwargs...
            ) where {C,T}

    ka = 0
    M = Helmholtz(sim)
    Σ0,Σ1,Σ2 = M.Σs
    f = M.sim.self_energy.f
    # f = map(f->(ω->f(ω)),M.sim.self_energy.f)
    Fs, fχs = _compute_Fsfχs(sim,1)
    As = vcat([Σ0,Σ1,Σ2,sim.laplacian.l0,M.αε],Fs)
    fs = [f[1],f[2],one,one,ω->ω^2,map(ϝ->(ω->ω^2*ϝ(ω)),fχs)...]
    nep = SPMF_NEP(As,fs;kwargs...)
    return HelmholtzNEP(M,nep,6:length(As))
end


"""
    MaxwellNEP(::Simulation{1}; ky=0, kz=0) -> nep
"""
function MaxwellNEP(
            sim::Simulation{1,C,T};
            ky::Number = 0,
            kz::Number = 0,
            kwargs...
            ) where {C,T}

    ka = 0
    M = Maxwell(sim; ky=ky, kz=kz)
    f = map(f->(ω->f(ω,ka,ky,kz)),M.sim.self_energy.f)
    Σ1,Σ2,Σ3 = M.Σs
    Fs, fχs = _compute_Fsfχs(sim,3)
    As = vcat([Σ1,Σ2,Σ3,sim.curlcurl(ky,kz),M.αε],Fs)
    fs = [f[1],f[2],one,one,ω->-ω^2,map(ϝ->(ω->-ω^2*ϝ(ω)),fχs)...]
    nep = SPMF_NEP(As,fs;kwargs...)
    return MaxwellNEP(M,nep,6:length(As))
end


################################################################################

function _compute_Fsfχs(sim::Simulation{1},m::Integer)
    χs = map(d->d.χ,sim.dispersive_domains)
    fχs = map(x->(ω->susceptability(x,ω)),χs)
    Fs = Vector{SparseMatrixCSC{ComplexF64,Int}}(undef,length(sim.dispersive_domains))
    for d ∈ eachindex(Fs) Fs[d] = spdiagm(0=>repeat(sim.Fs[d],m)) end
    return Fs,fχs
end

end # module

using .Spectral2D

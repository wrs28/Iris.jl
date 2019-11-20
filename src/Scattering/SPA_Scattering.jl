module SPA_Scattering

export scattering_spa

using ...Common
using LinearAlgebra
using PolynomialRoots
using NLsolve

function scattering_spa(sim::Simulation{1},ω::Number,b::Real)
    sct = scattering(sim,ω)
    cf = maxwell_cf(sim)
    cf(ω)
    φ = sct.tot
    ηs, us = maxwell_eigs(cf,ω; nev=5)
    _,ind = findmin(abs.(ηs))
    η, u = ηs[ind], us(ind)
    ηF = repeat(η*sim.F,3)
    DF = Vector(diag(cf.M.αχ))
    α = sum(u.^2 .*abs2.(u).*(ηF+DF))*sim.dx
    β = sum(u.^2 .*conj(u).*φ.*(ηF+2DF))*sim.dx
    γ = sum(u.^3 .*conj(φ).*(ηF+DF))*sim.dx
    δ = sum(u.^2 .*abs2.(φ).*(ηF+2DF))*sim.dx + η
    ε = sum(abs2.(u).*φ.^2 .*DF)*sim.dx
    ζ = sum(u.*φ.*abs2.(φ).*DF)*sim.dx
    as = [α,β,γ,δ,ε,ζ]
    pt = PhaseTuner(as,b,Vector{ComplexF64}(undef,3))
    results = nlsolve(pt,[0.])
    θ = mod2pi(results.zero[1])
end

struct PhaseTuner
    as::Vector{ComplexF64}
    b::Float64
    r::Vector{ComplexF64}
end
function (pt::PhaseTuner)(x)
    as = pt.as; b = pt.b
    α = as[1]; β = as[2]; γ = as[3]; δ = as[4]; ε = as[5]; ζ = as[6]
    b *= cis(x[1])
    a0 = ζ*b*abs2(b)
    a1 = δ*abs2(b)+ε*b^2+η
    a2 = β*conj(b)+γ*b
    a3 = α
    pt.r[:] = roots([a0,a1,a2,a3];polish=true)
    return pt.r
end
function (pt::PhaseTuner)(F,x)
    pt(x)
    _,ind = findmin(abs.(imag(pt.r)))
    F[1] = imag(pt.r[ind])
    return nothing
end


end # module

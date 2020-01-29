module TimeDomain1D

export _αβ

using ..Common
using ..DivGradCurls
import ..Common.DEFAULT_CFL_NUMBER
import LinearAlgebra: mul!

import ..HelmholtzWaveFields
import ..HelmholtzFDTD
import ..propagate!


# HelmholtzField constructor from simulation
"""
    HelmholtzWaveFields(simulation) -> hf
"""
HelmholtzWaveFields

function HelmholtzWaveFields(sim::Simulation{1,Common.Symmetric})
    φ = (ScalarField(sim.x, 1, Float64), ScalarField(sim.x, 1, Float64))
    ∇Φ = ScalarField(sim.x_half[1], 1, Float64)
    return HelmholtzWaveFields(φ, ((∇Φ,deepcopy(∇Φ)),))
end


"""
    HelmholtzFDTD(sim; dt=sim.dx*$DEFAULT_CFL_NUMBER)
"""
HelmholtzFDTD(sim::Simulation{1,Common.Symmetric}; dt::Real=sim.dx*DEFAULT_CFL_NUMBER) = HelmholtzFDTD(sim, dt)

function _αβ(sim::Simulation{1,Common.Symmetric}, dt::Real)
    σ = real(sim.σ[1])
    α = [(1 .- σ*dt/2)./(1 .+ σ*dt/2), dt ./(1 .+ σ*dt/2)]

    σ = real(sim.σ_half[1])
    β = [(1 .- σ*dt/2)./(1 .+ σ*dt/2), dt ./(1 .+ σ*dt/2)]
    return α, β
end

# @inline
function propagate!(fdtd::HelmholtzFDTD{1,2}, n::Integer=1)
    dt = fdtd.dt

    φ = fdtd.fields.φ
    ∇Φ = fdtd.fields.∇Φ

    grad = fdtd.grad
    div = fdtd.div

    α = fdtd.α
    β = fdtd.β

    foreach(1:n) do i
        # @fastmath

        mul!(φ[2], div, ∇Φ[1][1])
        # @fastmath @inbounds @simd
        for j ∈ eachindex(φ[1])
            φ[1][j] = α[1][j]*φ[1][j] + α[2][j]*φ[2][j]
        end
        # @fastmath @inbounds
        fdtd.t[1] += dt/2

        # @fastmath
        mul!(∇Φ[1][2], grad, φ[1])
        # @fastmath @inbounds @simd
        for j ∈ eachindex(∇Φ[1][1])
            ∇Φ[1][1][j] = β[1][j]*∇Φ[1][1][j] + β[2][j]*∇Φ[1][2][j]
        end
        # @fastmath @inbounds
        fdtd.t[1] += dt/2
    end
    fdtd.n[1] += n
    return nothing
end

end # module

using .TimeDomain1D

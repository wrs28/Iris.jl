module TimeDomain2D

export _αβ

using ..Common
using ..DivGradCurls
using Plots
using ProgressMeter
using RecipesBase

import ..Common.DEFAULT_CFL_NUMBER
import LinearAlgebra: mul!

import ..AbstractFDTD
import ..AbstractFields
import ..HelmholtzPointSource
import ..HelmholtzWaveFields
import ..HelmholtzFDTD
import ..propagate!


# HelmholtzField constructor from simulation
"""
    HelmholtzWaveFields(simulation) -> hf
"""
HelmholtzWaveFields

function HelmholtzWaveFields(sim::Simulation{2,Common.Symmetric})
    φ = (ScalarField(sim.x, 1, Float64), ScalarField(sim.x, 1, Float64))
    ∇Φₓ = ScalarField(sim.x_half[1], 1, Float64)
    ∇Φᵤ = ScalarField(sim.x_half[2], 1, Float64)
    return HelmholtzWaveFields(φ, ((∇Φₓ,deepcopy(∇Φₓ)),(∇Φᵤ,deepcopy(∇Φᵤ))))
end


"""
    HelmholtzFDTD(sim; source=0, dt=sim.dx*$DEFAULT_CFL_NUMBER, plotoptions...) -> fdtd
"""
function HelmholtzFDTD(
            sim::Simulation{2,Common.Symmetric};
            source = HelmholtzPointSource(sim, 0, 0, 0, 0),
            dt::Real = sim.dx*DEFAULT_CFL_NUMBER,
            kwargs...)

    return HelmholtzFDTD(sim, dt, source; kwargs...)
end

function _αβ(sim::Simulation{2,Common.Symmetric}, dt::Real)
    σ = real(sim.σ[1])
    α = [(1 .- σ*dt/2)./(1 .+ σ*dt/2), dt ./(1 .+ σ*dt/2)./sim.ε]

    σ = real(sim.σ_half[1])
    β = [(1 .- σ*dt/2)./(1 .+ σ*dt/2), dt ./(1 .+ σ*dt/2)]
    return α, β
end


"""
    HelmholtzPointSource(sim, xoft) -> ps
"""
function HelmholtzPointSource(sim::Simulation{2,Common.Symmetric}, xoft, ωoft, aoft=1, ϕoft=0) where F
    σ = 4sim.dx
    N = sqrt(2π)*σ
    return HelmholtzPointSource(xoft, ωoft, aoft, ϕoft, σ^2, N)
end

"""
    propagate!(fdtd, [n=1; animate=false, verbose=false])
"""
@inline function propagate!(fdtd::HelmholtzFDTD{2,2}, n::Integer=1; animate::Bool=false, verbose::Bool=false)
    dt = fdtd.dt

    φ = getfield(fdtd.fields,:φ)
    ∇Φ = fdtd.fields.∇Φ

    grad = fdtd.grad
    div = fdtd.div

    α = fdtd.α
    β = fdtd.β

    options = fdtd.options

    verbose ? pg = Progress(n) : nothing
    foreach(1:n) do i
        @fastmath mul!(φ[2], div, ∇Φ[1][1])
        @fastmath @inbounds @simd for j ∈ eachindex(φ[1])
            φ[1][j] = α[1][j]*φ[1][j] + α[2][j]*φ[2][j]
            iszero(fdtd.source.aoft(fdtd.t[1])) ? nothing : φ[1][j] -= dt*fdtd.source(φ[1].positions[j], fdtd.t[1])
        end
        @fastmath @inbounds fdtd.t[1] += dt/2

        @fastmath mul!(∇Φ[1][2], grad, φ[1])
        @fastmath @inbounds @simd for j ∈ eachindex(∇Φ[1][1])
            ∇Φ[1][1][j] = β[1][j]*∇Φ[1][1][j] + β[2][j]*∇Φ[1][2][j]
        end
        @fastmath @inbounds fdtd.t[1] += dt/2

        verbose ? next!(pg) : nothing

        if animate && iszero(mod(fdtd.n[1]+i-options.start,options.interval))
            if options.by ∈ (abs,abs2)
                options.ylims[1] = 0
                options.ylims[2] = max(options.ylims[2],maximum(options.by,φ[1]))
            else
                options.ylims[1] = min(options.ylims[1],minimum(options.by,φ[1]))
                options.ylims[2] = max(options.ylims[2],maximum(options.by,φ[1]))
                options.ylims[2] = max(abs(options.ylims[1]),abs(options.ylims[2]))
                options.ylims[1] = -options.ylims[2]
            end
            frame(options.animation, plot(fdtd, options.by; ylims=options.ylims, grid=false, options.plotoptions...))
        end
    end
    fdtd.n[1] += n
    return nothing
end


################################################################################
# Plotting

@recipe function f(field::HelmholtzWaveFields{2,2}, by::Function)
    @series begin field.φ, by end
end

@recipe function f(fdtd::HelmholtzFDTD{2,2}, by::Function)
    @series begin fdtd.fields, by end
end

end # module

using .TimeDomain2D

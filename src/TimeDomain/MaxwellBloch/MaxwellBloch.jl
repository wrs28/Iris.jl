"""
    module MaxwellBlochFDTD
"""
module MaxwellBlochFDTD

# export MaxwellBloch
# export MaxwellBlochField
export maxwell_bloch

using Formatting
using LinearAlgebra
using SparseArrays
using ...Common
using ..Curls
using ..PolarizationFields
using ..InversionFields
# using ..AuxilliaryFields

import ...Common.PRINTED_COLOR_DARK
import ...Common.PRINTED_COLOR_LIGHT
import ...Common.PRINTED_COLOR_NUMBER
import ...Common.TIME_TO_SPACE_STEP_RATIO

include("RealElectricFields.jl")
using .RealElectricFields

const RealMagneticField = RealElectricField

struct MaxwellBlochField{N}
    E::RealElectricField{N}
    B::RealMagneticField{N}
    P::PolarizationField{N}
    D::InversionField{N}
end

function MaxwellBlochField(sim::Simulation{1})
    E = RealElectricField(sim.x)
    posB = map(x->Point(x),Curls.expand_B_pos(sim))
    B = RealMagneticField(posB)
    P = PolarizationField(sim.x)
    D = InversionField(sim.x)
    MaxwellBlochField(E,B,P,D)
end

function Base.show(io::IO,mbf::MaxwellBlochField)
    n = length(mbf)
    printstyled(io,"MaxwellBlochField",color=PRINTED_COLOR_DARK)
    print(io," (")
    printstyled(io,n,color=PRINTED_COLOR_NUMBER)
    print(io," points)")
end

Base.length(mb::MaxwellBlochField) = length(mb.D)

struct MaxwellBloch{N,TS}
    fields::MaxwellBlochField{N}
    differentials::MaxwellBlochField{N}
    curl::Curl{N,Float64}
    n::Vector{Int}
    t::Vector{Float64}
    dt::Float64
    source::TS
    ρE::Vector{Float64}
    σEp1⁻¹::Vector{Float64}
    ρB::Vector{Float64}
    σBp1⁻¹::Vector{Float64}
    ε::Vector{Float64}
end

maxwell_bloch(args...; kwargs...) = MaxwellBloch(args...; kwargs...)

Base.length(mb::MaxwellBloch) = length(mb.fields)

function MaxwellBloch(sim::Simulation{1}, ω::Real; dt::Real=sim.dx*TIME_TO_SPACE_STEP_RATIO)
    fields = MaxwellBlochField(sim)
    differentials = MaxwellBlochField(sim)
    curl = Curl(sim)
    n = [0]
    t = [dt/2]
    sc = 0

    σ = real(sim.σ[1])
    σEp1⁻¹ = 1 ./(1 .+ σ*dt/2)
    ρE = (1 .- σ*dt/2).*σEp1⁻¹
    ρE = vcat(fill(0,length(ρE)),repeat(ρE,2))
    σEp1⁻¹ = vcat(fill(0,length(σEp1⁻¹)),repeat(σEp1⁻¹,2))

    σ = real(sim.σ_half[1])
    σBp1⁻¹ = 1 ./(1 .+ σ*dt/2)
    ρB = (1 .- σ*dt/2).*σBp1⁻¹
    ρB = vcat(fill(0,length(ρB)),repeat(ρB,2))
    σBp1⁻¹ = vcat(fill(0,length(σBp1⁻¹)),repeat(σBp1⁻¹,2))

    ε = repeat(sim.ε,3)
    return MaxwellBloch{1,typeof(sc)}(fields,differentials,curl,n,t,dt,sc,ρE,
                σEp1⁻¹,ρB,σBp1⁻¹,ε)
end

function Base.show(io::IO,mb::MaxwellBloch)
    n = length(mb)
    printstyled(io,"MaxwellBloch",color=PRINTED_COLOR_LIGHT)
    print(io," (")
    printstyled(io,n,color=PRINTED_COLOR_NUMBER)
    print(io," points, ")
    printstyled(io,"n=$(mb.n[1])",color=PRINTED_COLOR_NUMBER)
    print(io,", ")
    printstyled(io,"t=",fmt("1.2e",mb.t[1]),color=PRINTED_COLOR_NUMBER)
    print(io,")")
end

@inline function (mb::MaxwellBloch{1})(n::Int=1)
    dt = mb.dt

    E = mb.fields.E
    dE = mb.differentials.E
    curlE = mb.curl.curlE

    B = mb.fields.B
    dB = mb.differentials.B
    curlB = mb.curl.curlB

    foreach(1:n) do i
        @fastmath mul!(dE.val,curlB,B.val)
        @fastmath @inbounds @simd for j ∈ eachindex(dE.val) E[j] = mb.ρE[j]*E[j] - mb.ε[j]\mb.σEp1⁻¹[j]*dt*dE[j] end
        @fastmath @inbounds mb.t[1] += dt/2

        @fastmath mul!(dB.val,curlE,E.val)
        @fastmath @inbounds @simd for j ∈ eachindex(dB.val) B[j] = mb.ρB[j]*B[j] + mb.σBp1⁻¹[j]*dt*dB[j] end
        @fastmath @inbounds mb.t[1] += dt/2
    end
    mb.n[1] += n
    return nothing
end

end # module

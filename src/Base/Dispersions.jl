module Dispersions

export AbstractDispersion,
TwoLevelSystem

abstract type AbstractDispersion end

"""
    tls = TwoLevelSystem(tls; :key1 => value1, :key2 => value2, ...)

new tls object from old, with modified fields
"""
struct TwoLevelSystem <: AbstractDispersion
    D₀::Float64
    k₀::Float64
    γp::Float64
    ω₀::Float64

    TwoLevelSystem(D₀=0, k₀=Inf, γp=1e8, ω₀=k₀) = new(D₀, k₀, γp, ω₀)
    TwoLevelSystem(tls::TwoLevelSystem; D₀=tls.D₀, k₀=tls.k₀, γp=tls.γp, ω₀=tls.ω₀) = TwoLevelSystem(tls.D₀, tls.k₀, tls.γp, tls.ω₀)

    function Base.show(io::IO, tls::TwoLevelSystem)
        if !get(io, :sub, false)
            print(io, "Two Level System:\n")
        end
        print(io, "\tD₀: ", tls.D₀, "\n",
        "\tω₀: ", tls.ω₀, "\n",
        "\tγ⟂: ", tls.γp)
    end
end

(tls::TwoLevelSystem)(ω::Number) = tls.γp/(ω-tls.ω₀+1im*tls.γp)
function (tls::TwoLevelSystem)(ω::Number,ωs::Array,ψ::Array)
    γω = tls(ω)
    h = zeros(Float64,size(ψ,1))
    for i ∈ eachindex(ωs)
        h += abs2(ωs[i])*abs2.(ψ[:,i])
    end
    return γω./(1 .+ h)
end

end # module

#NOTE: \perp does not work, reserved as infix operator
# use \bigbot instead

"""
    tls = TwoLevelSystem(tls; :key1 => value1, :key2 => value2, ...)

new tls object from old, with modified fields
"""
struct TwoLevelSystem <: AbstractDispersion
    ωₐ::Float64
    D₀::Float64
    γ⟘::Float64

    TwoLevelSystem(ωₐ, D₀ = 0, γ⟘ = 1e8) = new(ωₐ, D₀, γ⟘)
end
TwoLevelSystem(tls::TwoLevelSystem; ωₐ=tls.ωₐ, D₀=tls.D₀, γ=tls.γ⟘) = TwoLevelSystem(ωₐ, D₀, γ)

function Base.show(io::IO, tls::TwoLevelSystem)
    if !get(io, :sub, false)
        print(io, "Two Level System:\n")
    end
    print(io, "\tωₐ: ", tls.ωₐ, "\n",
    "\tD₀: ", tls.D₀, "\n",
    "\tγ⟘: ", tls.γ⟘)
end

(tls::TwoLevelSystem)(ω::Number) = tls.γ⟘/(ω-tls.ωₐ+1im*tls.γ⟘)
function (tls::TwoLevelSystem)(ω::Number,ωs::Array,ψ::Array)
    h = zeros(Float64,size(ψ,1))
    for i ∈ eachindex(ωs)
        h += abs2.(tls(ωs[i])*ψ[:,i])
    end
    return tls(ω)./(1 .+ h)
end

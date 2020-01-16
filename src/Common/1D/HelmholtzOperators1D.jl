function Helmholtz(
            sim::TSIM;
            m::Int = 1) where TSIM<:Simulation{1,Symmetric}

    ka, kb, kc = 0, 0, 0
    laplacian = sim.laplacian.l0

    Σs = _compute_Σs(sim)
    f = sim.self_energy.f
    D² = laplacian - Σs[1]*f[1](0,ka) - Σs[2]*f[2](0,ka) - Σs[3]

    αε = sim.laplacian.complex_scaling*spdiagm(0=>sim.ε)
    αεpFχ = copy(αε)
    Fχ = spdiagm(0=>zeros(ComplexF64,size(αε,1)))

    Fχs = Vector{SparseMatrixCSC{ComplexF64,Int}}(undef,m)
    dFχdψr = zeros(ComplexF64,length(sim),m,m)
    dFχdψi = copy(dFχdψr)
    dFχdω = copy(dFχdψi)
    dFχdϕ = copy(dFχdω)
    for μ ∈ eachindex(Fχs) Fχs[μ] = spdiagm(0=>zeros(ComplexF64,size(αε,1))) end
    Helmholtz{1,Symmetric,typeof(Σs),TSIM}(D²+I,D²,αεpFχ,laplacian,Σs,αε,Fχ,Fχs,dFχdψr,dFχdψi,dFχdω,dFχdϕ,sim,ka,kb,kc)
end


@inline function (h::Helmholtz{1})(ω::Number,args...)
    f = h.sim.self_energy.f
    ka = 0
    fL = f[1](ω,ka)
    fR = f[2](ω,ka)

    rows = rowvals(h.D²)
    vals = nonzeros(h.D²)
    _, n = size(h.D²)
    @inbounds for i ∈ 1:n
        col = i
        @fastmath @inbounds @simd for j ∈ nzrange(h.D², i)
            row = rows[j]
            vals[j] = h.laplacian[row,col] - h.Σs[1][row,col]*fL - h.Σs[2][row,col]*fR - h.Σs[3][row,col]
        end
    end

    helmholtz_susceptability!(h,ω,args...)
    ω² = ω^2

    rows = rowvals(h.A)
    vals = nonzeros(h.A)
    @inbounds for col ∈ 1:n
        @fastmath @inbounds @simd for j ∈ nzrange(h.A, col)
            row = rows[j]
            vals[j] = h.D²[row,col] + ω²*h.αεpFχ[row,col]
        end
    end
    return h.A
end

################################################################################

function _compute_Σs(sim::Simulation{1})
    Σ0L = sim.self_energy.Σ0[1]
    Σ0R = sim.self_energy.Σ0[2]
    # Σ1L = sim.self_energy.Σ1[1]
    # Σ1R = sim.self_energy.Σ1[2]
    return (Σ0L, Σ0R, sim.self_energy.Σ2)
end

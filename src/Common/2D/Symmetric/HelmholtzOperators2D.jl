function Helmholtz(
            sim::TSIM;
            m::Int = 1) where TSIM<:Simulation{2}

    ka, kb, kc = 0, 0, 0, 0
    laplacian = sim.laplacian.l0

    Σ01,Σ02,Σ11,Σ12,Σ2 = _compute_Σs(sim)
    Σs = (Σ01 + Σ11, Σ02 + Σ12, Σ2)
    f = sim.self_energy.f
    D² = Σs[1]*f[1](0,ka) + Σs[2]*f[2](0,ka) + Σs[3] + laplacian

    αε = _compute_αε(sim)
    αεpFχ = copy(αε)
    Fχ = spdiagm(0=>zeros(ComplexF64,size(αε,1)))

    Fχs = Vector{SparseMatrixCSC{ComplexF64,Int}}(undef,m)
    dFχdψr = zeros(ComplexF64,length(sim),m,m)
    dFχdψi = copy(dFχdψr)
    dFχdω = copy(dFχdψi)
    dFχdϕ = copy(dFχdω)
    for μ ∈ eachindex(Fχs) Fχs[μ] = spdiagm(0=>zeros(ComplexF64,size(αε,1))) end
    Helmholtz{1,typeof(Σs),TSIM}(D²+I,D²,αεpFχ,laplacian,Σs,αε,Fχ,Fχs,dFχdψr,dFχdψi,dFχdω,dFχdϕ,sim,ka,kb,kc)
end


@inline function (m::Helmholtz{1})(ω::Number,args...)
    f = m.sim.self_energy.f
    ka = 0
    f1 = f[1](ω,ka)
    f2 = f[2](ω,ka)

    rows = rowvals(m.D²)
    vals = nonzeros(m.D²)
    _, n = size(m.D²)
    @inbounds for i ∈ 1:n
        col = i
        @fastmath @inbounds @simd for j ∈ nzrange(m.D², i)
            row = rows[j]
            vals[j] = m.Σs[1][row,col]*f1 + m.Σs[2][row,col]*f2 + m.Σs[3][row,col] + m.laplacian[row,col]
        end
    end

    helmholtz_susceptability!(m,ω,args...)
    ω² = ω^2

    rows = rowvals(m.A)
    vals = nonzeros(m.A)
    @inbounds for col ∈ 1:n
        @fastmath @inbounds @simd for j ∈ nzrange(m.A, col)
            row = rows[j]
            vals[j] = m.D²[row,col] - ω²*m.αεpFχ[row,col]
        end
    end
    return m.A
end



################################################################################

# linear
@inline function helmholtz_susceptability!(m::Helmholtz{1},ω::Number)
    sim = m.sim
    @fastmath @inbounds @simd for i ∈ eachindex(m.Fχ.nzval)
        j = mod1(i,length(sim))
        χ = susceptability(sim.χ[j], ω)
        m.Fχ.nzval[i] = sim.F[j]*χ
    end
    @fastmath @inbounds @simd for i ∈ eachindex(m.αε.nzval) m.αεpFχ.nzval[i] = m.αε.nzval[i] + m.Fχ.nzval[i] end
    return nothing
end


# nonlinear + local jacobian
@inline function helmholtz_susceptability!(m::Helmholtz{1},ω,ωs::Vector,ψs::ElectricField)
    sim = m.sim
    foreach(d->susceptability(d.χ,ω,ωs,ψs), sim.dispersive_domains)
    @inbounds for μ ∈ eachindex(m.Fχs)
        @inbounds for i ∈ eachindex(m.Fχ.nzval)
            j = mod1(i,length(sim))
            if typeof(sim.χ[j])<:TwoLevelSystem
                χ::TwoLevelSystem = sim.χ[j]
                μ==1 ? m.Fχ.nzval[i] = sim.F[j]*χ.chi[i] : nothing
                m.Fχs[μ].nzval[i] = sim.F[j]*χ.chis[i,μ]
                @inbounds @simd for ν ∈ eachindex(m.Fχs)
                    m.dFχdψr[i,μ,ν] = sim.F[j]*χ.dχdψr[i,μ,ν]
                    m.dFχdψi[i,μ,ν] = sim.F[j]*χ.dχdψi[i,μ,ν]
                    m.dFχdω[i,μ,ν] = sim.F[j]*χ.dχdω[i,μ,ν]
                    m.dFχdϕ[i,μ,ν] = sim.F[j]*χ.dχdϕ[i,μ,ν]
                end
            end
        end
    end
    @fastmath @inbounds @simd for i ∈ eachindex(m.αε.nzval) m.αεpFχ.nzval[i] = m.αε.nzval[i] + m.Fχ.nzval[i] end
    return nothing
end


################################################################################

function _compute_Σs(sim::Simulation{1})
    Σ01 = sim.self_energy.Σ0[1]
    Σ02 = sim.self_energy.Σ0[2]
    Σ11 = sim.self_energy.Σ1[1]
    Σ12 = sim.self_energy.Σ1[2]
    return Σ01,Σ02,Σ11,Σ12,sim.self_energy.Σ2
end


_compute_αε(sim::Simulation{1}) = spdiagm(0=>sim.α[1].*sim.ε)

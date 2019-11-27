function Maxwell(
            sim::TSIM;
            ky::Number=0,
            kz::Number=0,
            m::Int = 1) where TSIM<:Simulation{1}

    kx, ka, kb, kc = 0, 0, 0, 0
    curlcurl = sim.curlcurl(ky,kz)

    Σ01,Σ02,Σ11,Σ12,Σ2 = _compute_Σs(sim,ky,kz)
    Σs = (Σ01 + Σ11, Σ02 + Σ12, Σ2)
    f = sim.self_energy.f
    D² = Σs[1]*f[1](0,ka,ky,kz) + Σs[2]*f[2](0,ka,ky,kz) + Σs[3] + curlcurl

    αε = _compute_αε(sim)
    αεpFχ = copy(αε)
    Fχ = spdiagm(0=>zeros(ComplexF64,size(αε,1)))

    Fχs = Vector{SparseMatrixCSC{ComplexF64,Int}}(undef,m)
    dFχdψr = zeros(ComplexF64,3length(sim),m,m)
    dFχdψi = copy(dFχdψr)
    dFχdω = copy(dFχdψi)
    dFχdϕ = copy(dFχdω)
    for μ ∈ eachindex(Fχs) Fχs[μ] = spdiagm(0=>zeros(ComplexF64,size(αε,1))) end
    Maxwell{1,typeof(Σs),TSIM}(D²+I,D²,αεpFχ,curlcurl,Σs,αε,Fχ,Fχs,dFχdψr,dFχdψi,dFχdω,dFχdϕ,sim,kx,ky,kz,ka,kb,kc)
end


@inline function (m::Maxwell{1})(ω::Number,args...)
    ka, ky, kz = m.ka, m.ky, m.kz
    f = m.sim.self_energy.f

    f1 = f[1](ω,ka,ky,kz)
    f2 = f[2](ω,ka,ky,kz)

    rows = rowvals(m.D²)
    vals = nonzeros(m.D²)
    _, n = size(m.D²)
    @inbounds for i ∈ 1:n
        col = i
        @fastmath @inbounds @simd for j ∈ nzrange(m.D², i)
            row = rows[j]
            vals[j] = m.Σs[1][row,col]*f1 + m.Σs[2][row,col]*f2 + m.Σs[3][row,col] + m.curlcurl[row,col]
        end
    end

    maxwell_susceptability!(m,ω,args...)
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
@inline function maxwell_susceptability!(m::Maxwell{1},ω::Number)
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
@inline function maxwell_susceptability!(m::Maxwell{1},ω,ωs::Vector,ψs::ElectricField)
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

@inline function _compute_Σs(sim::Simulation{1},ky::Number,kz::Number)

    ky², kz², kykz = ky^2, kz^2, ky*kz

    Σ0 = sim.self_energy.Σ0
    Σ1 = sim.self_energy.Σ1
    Σ2 = sim.self_energy.Σ2

    I1 = sparse([2,3],[2,3],[1,1],3,3)
    I2 = sparse([2,3,1,1],[1,1,2,3],[ky,kz,ky,kz],3,3)

    Σ01 = kron(I1,Σ0[1])
    Σ02 = kron(I1,Σ0[2])
    Σ11 = kron(I2,Σ1[1])
    Σ12 = kron(I2,Σ1[2])
    Σ2 = kron(sparse([1,2,3,2,3],[1,2,2,3,3],[ky²+kz²,kz²,kykz,kykz,ky²],3,3),Σ2)
    return Σ01,Σ02,Σ11,Σ12,Σ2
end


_compute_αε(sim::Simulation{1})::SparseMatrixCSC{ComplexF64,Int} =
    kron(I3,spdiagm(0=>sim.α[1].*sim.ε))

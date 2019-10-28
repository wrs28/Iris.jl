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
    αεpχ = similar(αε)
    αχ = similar(αε)

    αχs = Vector{SparseMatrixCSC{ComplexF64,Int}}(undef,m)
    αdχdψr = similar(αχs)
    αdχdψi = similar(αdχdψr)
    for i ∈ eachindex(αχs) αχs[i] = similar(αε) end
    for i ∈ eachindex(αdχdψr) αdχdψr[i] = similar(αε) end
    for i ∈ eachindex(αdχdψr) αdχdψi[i] = similar(αε) end

    Maxwell{1,typeof(Σs),TSIM}(D²+αχ,D²,αεpχ,curlcurl,Σs,αε,αχ,αχs,αdχdψr,αdχdψi,sim,kx,ky,kz,ka,kb,kc)
end


function (m::Maxwell{1})(ω::Number,args...)
    ka, ky, kz = m.ka, m.ky, m.kz
    f = m.sim.self_energy.f

    f1 = f[1](ω,ka,ky,kz)
    f2 = f[2](ω,ka,ky,kz)

    rows = rowvals(m.D²)
    vals = nonzeros(m.D²)
    _, n = size(m.D²)
    for i ∈ 1:n
        col = i
        for j ∈ nzrange(m.D², i)
            row = rows[j]
            vals[j] = m.Σs[1][row,col]*f1 + m.Σs[2][row,col]*f2 + m.Σs[3][row,col] + m.curlcurl[row,col]
        end
    end

    maxwell_susceptability!(m,ω,args...)
    ω² = ω^2

    rows = rowvals(m.A)
    vals = nonzeros(m.A)
    _, n = size(m.A)
    for i ∈ 1:n
        col = i
        for j ∈ nzrange(m.A, i)
            row = rows[j]
            vals[j] = m.D²[row,col] - ω²*m.αεpχ[row,col]
        end
    end

    return m.A
end



################################################################################

# linear
function maxwell_susceptability!(m::Maxwell{1},ω::Number,args...)
    sim = m.sim
    n::Int = length(sim)
    χs = map(d->d.χ,sim.domains)
    Fs = map(d->d.pump,sim.domains)
    χχ = map(x->susceptability(x,ω,args...),χs)
    foreach((χ,d)->begin
                    for i ∈ eachindex(m.αχ.nzval)
                        j = mod1(i,n)
                        if d==m.sim.domain_indices[j]
                            m.αχ.nzval[i]=0
                            for x ∈ χ m.αχ.nzval[i] += Fs[d](m.sim.x[j])*m.sim.α[1][j]*x end
                        end
                    end
                end, χχ, eachindex(χχ))
    for i ∈ eachindex(m.αε.nzval) m.αεpχ.nzval[i] = m.αε.nzval[i] + m.αχ.nzval[i] end
    return nothing
end


# nonlinear + local jacobian
function maxwell_susceptability!(m::Maxwell{1},ωs::Vector,ψs::ElectricField)
    sim = m.sim
    n::Int = length(sim)
    χs = map(d->d.χ,sim.domains)
    Fs = map(d->d.pump,sim.domains)
    foreach(x->susceptability(x,ωs,ψs),χs)
    foreach((χ,d)->begin
                    for μ ∈ eachindex(ml.αχ)
                        for i ∈ eachindex(ml.αχ[μ].nzval)
                            j = mod1(i,n)
                            if d==ml.sim.domain_indices[j]
                                ml.αχ[μ].nzval[i]=0
                                for x ∈ χ
                                    ml.αχ[μ].nzval[i] += Fs[d](ml.sim.x[j])*ml.sim.α[1][j]*x.chi[i,μ]
                                    ml.αdχdψr[μ].nzval[i] += Fs[d](ml.sim.x[j])*ml.sim.α[1][j]*x.dχdψr[i,μ]
                                    ml.αdχdψi[μ].nzval[i] += Fs[d](ml.sim.x[j])*ml.sim.α[1][j]*x.dχdψi[i,μ]
                                end
                            end
                        end
                    end
                end, χs, eachindex(χs))

    for μ ∈ eachindex(ml.αεχ) for i ∈ eachindex(ml.αε.nzval) ml.αεχ[μ].nzval[i] = ml.αε.nzval[i] + ml.αχ[μ].nzval[i] end end
    return nothing
end


################################################################################

function _compute_Σs(sim::Simulation{1},ky::Number,kz::Number)

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

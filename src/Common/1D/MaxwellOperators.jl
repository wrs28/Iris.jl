
maxwell_lep(sim::Simulation{1};ky::Number=0,kz::Number=0) = Maxwell_LEP(sim;ky=ky,kz=kz)

function Maxwell_LEP(
            sim::Simulation{1};
            ky::Number=0,
            kz::Number=0)

    αε = _compute_αε(sim)
    return Maxwell_LEP(sim,
            _compute_Σs(sim,ky,kz),
            sim.curlcurl(ky,kz),
            αε,
            spdiagm(0=>zeros(ComplexF64,size(αε,1))),
            0,ky,kz,0,0,0)
end
function (ml::Maxwell_LEP{1})(ω::Number)
    sim = ml.sim
    n::Int = length(sim.x)
    ka, ky, kz = ml.ka, ml.ky, ml.kz
    f = ml.sim.self_energy.f
    Σ01,Σ02,Σ11,Σ12,Σ2 = ml.Σs

    Σ0 = Σ01*f[1](ω,ka,ky,kz) + Σ02*f[2](ω,ka,ky,kz)
    Σ1 = Σ11*f[1](ω,ka,ky,kz) + Σ12*f[2](ω,ka,ky,kz)

    χs = map(d->d.χ,sim.domains)
    Fs = map(d->d.pump,sim.domains)
    χχ = map(x->susceptability(x,ω),χs)
    foreach((χ,d)->begin
                    for i ∈ eachindex(ml.αχ.nzval)
                        j = mod1(i,n)
                        if d==ml.sim.domain_indices[j]
                            ml.αχ.nzval[i]=0
                            for x ∈ χ
                                ml.αχ.nzval[i] += Fs[d](ml.sim.x[j])*ml.sim.α[1][j]*x
                            end
                        end
                    end
                end, χχ, eachindex(χχ))
    return Σ0 + Σ1 + Σ2 + ml.curlcurl, ml.αε + ml.αχ
end
function (ml::Maxwell_LEP{1})(ω::Number,ωs::Array,ψs::Array)
    sim = ml.sim
    n::Int = length(sim.x)
    ka, ky, kz = ml.ka, ml.ky, ml.kz
    f = ml.sim.self_energy.f
    Σ01,Σ02,Σ11,Σ12,Σ2 = ml.Σs

    Σ0 = Σ01*f[1](ω,ka,ky,kz) + Σ02*f[2](ω,ka,ky,kz)
    Σ1 = Σ11*f[1](ω,ka,ky,kz) + Σ12*f[2](ω,ka,ky,kz)

    χs = map(d->d.χ,sim.domains)
    Fs = map(d->d.pump,sim.domains)
    foreach(x->susceptability(x,ω,ωs,ψs),χs)
    foreach((χ,d)->begin
                    for i ∈ eachindex(ml.αχ.nzval)
                        j = mod1(i,n)
                        if d==ml.sim.domain_indices[j]
                            ml.αχ.nzval[i]=0
                            for x ∈ χ
                                ml.αχ.nzval[i] += Fs[d](ml.sim.x[j])*ml.sim.α[1][j]*x.chi[i]
                            end
                        end
                    end
                end, χs, eachindex(χs))
    return Σ0 + Σ1 + Σ2 + ml.curlcurl, ml.αε + ml.αχ
end

maxwell(sim::Simulation{1};ky::Number=0,kz::Number=0) = Maxwell(maxwell_lep(sim))
function (m::Maxwell{1})(ω::Number,args...)
    A,B =m.lep(ω,args...)
    return A-B*ω^2
end

function maxwell_nep(
            sim::Simulation{1,C,T};
            ky::Number = 0,
            kz::Number = 0,
            FType::DataType = T,
            kwargs...
            ) where {C,T}

    ka = 0
    f = map(f->(ω->f(ω,ka,ky,kz)),sim.self_energy.f)
    Σ01,Σ02,Σ11,Σ12,Σ2 = _compute_Σs(sim,ky,kz)
    αε = _compute_αε(sim)
    αF, fχ = _compute_αχ(sim)
    As = vcat([Σ01,Σ02,Σ11,Σ12,Σ2,sim.curlcurl(ky,kz),αε],αF)
    fs = [f[1],f[2],f[1],f[2],one,one,ω->-ω^2,map(ϝ->(ω->-ω^2*ϝ(ω)),fχ)...]
    return SPMF_NEP(As,fs;kwargs...)
end


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


function _compute_αχ(sim::Simulation{1})
    χs = flatten(map(d->d.χ,sim.domains),AbstractDispersion)
    domain_inds = flatten(map(d->ntuple(identity,length(d.χ)),sim.domains))
    fχs = map(x->(ω->susceptability(x,ω)),χs)

    n = length(sim.x)
    Fs = map(d->sim.domains[d].pump,domain_inds)
    αFs = Array{SparseMatrixCSC{ComplexF64,Int},1}(undef,length(domain_inds))
    for d ∈ eachindex(αFs)
        αFs[d] = spdiagm(0=>zeros(ComplexF64,3n))
        for i ∈ 1:3n
            j = mod1(i,n)
            if domain_inds[d]==sim.domain_indices[j]
                αFs[d].nzval[i] = Fs[d](sim.x[j])*sim.α[1][j]
            end
        end
    end
    return αFs,fχs
end

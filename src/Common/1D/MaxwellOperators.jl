function maxwell_lep(sim::Simulation{1};ky::Number=0,kz::Number=0)

    ka = 0

    Σ01,Σ02,Σ11,Σ12,Σ2 = _compute_Σs(sim,ky,kz)
    αε = _compute_αε(sim)
    f = sim.self_energy.f

    f0(ω) = Σ01*f[1](ω,ka,ky,kz) + Σ02*f[2](ω,ka,ky,kz)
    f1(ω) = Σ11*f[1](ω,ka,ky,kz) + Σ12*f[2](ω,ka,ky,kz)
    f2(ω) = Σ2

    return ω -> (f0(ω) + f1(ω) + f2(ω) + sim.curlcurl(ky,kz), -αε)
end

function maxwell(sim::Simulation{1};ky::Number=0,kz::Number=0)
    return ω -> begin
                A,B = maxwell_lep(sim;ky=ky,kz=kz)(ω)
                return A-B*ω^2
            end
end


function maxwell_nep(
            sim::Simulation{1,C,T};
            ky::Number = 0,
            kz::Number = 0,
            FType::DataType = T,
            kwargs...
            ) where {C,T}

    ka = 0
    Σ01,Σ02,Σ11,Σ12,Σ2 = _compute_Σs(sim,ky,kz)
    f = map(f->(ω->f(ω,ka,ky,kz)),sim.self_energy.f)
    αε = _compute_αε(sim)
    As = [Σ01,Σ02,Σ11,Σ12,Σ2,sim.curlcurl(ky,kz),αε]
    fs = [f[1],f[2],f[1],f[2],one,one,ω->ω^2]
    return SPMF_NEP(As,fs;kwargs...)
end


function _compute_Σs(sim::Simulation{1},ky::Number,kz::Number)

    ky² = ky^2
    kz² = kz^2
    kykz = ky*kz

    Σ0 = sim.self_energy.Σ0
    Σ1 = sim.self_energy.Σ1
    Σ2 = sim.self_energy.Σ2

    Σ01 = kron(sparse([2,3],[2,3],[1,1],3,3),Σ0[1])
    Σ02 = kron(sparse([2,3],[2,3],[1,1],3,3),Σ0[2])
    Σ11 = kron(sparse([2,3,1,1],[1,1,2,3],[ky,kz,ky,kz],3,3),Σ1[1])
    Σ12 = kron(sparse([2,3,1,1],[1,1,2,3],[ky,kz,ky,kz],3,3),Σ1[2])
    Σ2 = kron(sparse([1,2,3,2,3],[1,2,2,3,3],[ky²+kz²,kz²,kykz,kykz,ky²],3,3),Σ2)
    return Σ01,Σ02,Σ11,Σ12,Σ2
end

_compute_αε(sim::Simulation{1})::SparseMatrixCSC{ComplexF64,Int} =
    kron(sparse(I,3,3),spdiagm(0=>sim.α[1].*sim.ε))

module SelfEnergy1D

using ..BoundaryConditions
using ..Domains
using ..Points
using SparseArrays

import ..Symmetric, ..Unsymmetric
import LinearAlgebra: I
import ..SelfEnergy

function SelfEnergy{Symmetric}(domain::LatticeDomain{1,Symmetric}, α_half)

    N = length(α_half)-1
    a = domain.shape.a
    α_half⁻¹ = 1 ./α_half
    dx = domain.dx

    bc1 = domain.boundary.bcs[1]
    if typeof(bc1)<:DirichletBC{1}
        Σ0L = sparse([1],[1],[-α_half⁻¹[1]/dx^2],N,N)
        Σ1L = sparse([1],[1],[+1/dx],N,N)
        fL = FΣ1D(dx,0)
    elseif typeof(bc1)<:NeumannBC{1}
        Σ0L = sparse([1],[1],[+α_half⁻¹[1]/dx^2],N,N)
        Σ1L = sparse([1],[1],[-1/dx],N,N)
        fL = FΣ1D(dx,0)
    elseif typeof(bc1)<:FloquetBC{1}
        # Σ0L = sparse([li],[ri],[α_half⁻¹[ri+1]/dx^2],N,N)
        # Σ1L = sparse([li],[ri],[-1/dx],N,N)
        # fL = (ω,ka,_...) -> ffl(ka,a)
        throw("not yet")
    elseif typeof(bc1)<:MatchedBC{1}
        Σ0L = sparse([1],[1],[+α_half⁻¹[1]/dx^2],N,N)
        Σ1L = sparse([1],[1],[-1/dx],N,N)
        if isempty(bc1.in)
            fL = FΣ1D(dx,+1)
        else
            fL = FΣ1D(dx,-1)
        end
    else
        Σ0L = spzeros(N,N)
        Σ1L = spzeros(N,N)
        fL = FΣ1D(dx,0)
    end

    bc2 = domain.boundary.bcs[2]
    if typeof(bc2)<:DirichletBC{2}
        Σ0R = sparse([N],[N],[-α_half⁻¹[N+1]/dx^2],N,N)
        Σ1R = sparse([N],[N],[-1/dx],N,N)
        fR = FΣ1D(dx,0)
    elseif typeof(bc2)<:NeumannBC{2}
        Σ0R = sparse([N],[N],[+α_half⁻¹[N+1]/dx^2],N,N)
        Σ1R = sparse([N],[N],[+1/dx],N,N)
        fR = FΣ1D(dx,0)
    elseif typeof(bc2)<:FloquetBC{2}
        # Σ0R = sparse([ri],[1],[-α_half⁻¹[li]/dx^2],N,N)
        # Σ1R = sparse([ri],[1],[+1/dx],N,N)
        # fR = (ω,ka,_...) -> ffr(ka,a)
        throw("not yet")
    elseif typeof(bc2)<:MatchedBC{2}
        Σ0R = sparse([N],[N],[+α_half⁻¹[N+1]/dx^2],N,N)
        Σ1R = sparse([N],[N],[+1/dx],N,N)
        if isempty(bc2.in)
            fR = FΣ1D(dx,+1)
        else
            fR = FΣ1D(dx,-1)
        end
    else
        Σ0R = spzeros(N,N)
        Σ1R = spzeros(N,N)
        fR = FΣ1D(dx,0)
    end
    Σ2 = spzeros(N,N)

    fs = (fL,fR)

    return SelfEnergy{1,Symmetric,2,typeof(fs)}((-Σ0L,-Σ0R),(Σ1L,Σ1R),Σ2,(fL,fR))
end

struct FΣ1D <: Function
    dx::Float64
    sign::Int
end
function (fm::FΣ1D)(ω::Number,ka=0,ky=0,kz=0)
    if iszero(fm.sign)
        return complex(1.0,0)
    else
        k̂ = sqrt(ω^2 - ky^2 - kz^2)*fm.dx/2
        abs(k̂)<1 || throw(DomainError(fm.dx,"lattice spacing $(fm.dx) too large. for |ω|=$(abs(ω)), ky=$ky, kz=$kz, should be < $(abs(1/sqrt(ω^2 - ky^2 - kz^2)))"))
        kx = 2asin(k̂)/fm.dx
        return exp(fm.sign*1im*kx*fm.dx)
    end
end
function (fm::FΣ1D)(ω::Matrix,ka=0,ky=0,kz=0)
    if iszero(fm.sign)
        return one(convert(Matrix{ComplexF64},ω))
    else
        kx::Matrix{ComplexF64} = 2asin(sqrt(ω^2 - one(ω)*ky^2 - one(ω)*kz^2)*fm.dx/2)/fm.dx
        return exp(fm.sign*1im*kx*fm.dx)
    end
end

end

using .SelfEnergy1D

module SelfEnergy2DSymmetric

using ..BoundaryConditions
using ..Domains
using ..Points
using SparseArrays

import ..Symmetric, ..Unsymmetric
import LinearAlgebra: I
import ..SelfEnergy

function SelfEnergy{Symmetric}(domain::LatticeDomain{2,Symmetric,Cartesian}, αx_half::Vector, αy_half::Vector)

    Nx = length(αx_half)-1
    Ny = length(αy_half)-1
    N = Nx*Ny

    Ix = sparse(I,Nx,Nx)
    Iy = sparse(I,Ny,Ny)
    a = domain.shape.a
    αx_half⁻¹ = 1 ./αx_half
    αy_half⁻¹ = 1 ./αy_half
    dx = domain.dx
    dy = domain.dy

    bc1 = domain.boundary.bcs[1]
    if typeof(bc1)<:DirichletBC{1}
        Σ0L = kron(Iy,sparse([1],[1],[-αx_half⁻¹[1]/dx^2],Nx,Nx))
        Σ1L = kron(Iy,sparse([1],[1],[+1/dx],Nx,Nx))
        fL = FΣ2D(dx,0)
    elseif typeof(bc1)<:NeumannBC{1}
        Σ0L = kron(Iy,sparse([1],[1],[+αx_half⁻¹[1]/dx^2],Nx,Nx))
        Σ1L = kron(Iy,sparse([1],[1],[-1/dx],Nx,Nx))
        fL = FΣ2D(dx,0)
    elseif typeof(bc1)<:FloquetBC{1}
        # Σ0L = sparse([li],[ri],[α_half⁻¹[ri+1]/dx^2],N,N)
        # Σ1L = sparse([li],[ri],[-1/dx],N,N)
        # fL = (ω,ka,_...) -> ffl(ka,a)
        throw("not yet")
    elseif typeof(bc1)<:MatchedBC{1}
        Σ0L = kron(Iy,sparse([1],[1],[+αx_half⁻¹[1]/dx^2],Nx,Nx))
        Σ1L = kron(Iy,sparse([1],[1],[-1/dx],Nx,Nx))
        if isempty(bc1.in)
            fL = FΣ2D(dx,+1)
        else
            fL = FΣ2D(dx,-1)
        end
    else
        Σ0L = spzeros(N,N)
        Σ1L = spzeros(N,N)
        fL = FΣ2D(dx,0)
    end

    bc2 = domain.boundary.bcs[2]
    if typeof(bc2)<:DirichletBC{2}
        Σ0R = kron(Iy,sparse([Nx],[Nx],[-αx_half⁻¹[Nx+1]/dx^2],Nx,Nx))
        Σ1R = kron(Iy,sparse([Nx],[Nx],[-1/dx],Nx,Nx))
        fR = FΣ2D(dx,0)
    elseif typeof(bc2)<:NeumannBC{2}
        Σ0R = kron(Iy,sparse([Nx],[Nx],[+αx_half⁻¹[Nx+1]/dx^2],Nx,Nx))
        Σ1R = kron(Iy,sparse([Nx],[Nx],[+1/dx],Nx,Nx))
        fR = FΣ2D(dx,0)
    elseif typeof(bc2)<:FloquetBC{2}
        # Σ0R = sparse([ri],[1],[-α_half⁻¹[li]/dx^2],N,N)
        # Σ1R = sparse([ri],[1],[+1/dx],N,N)
        # fR = (ω,ka,_...) -> ffr(ka,a)
        throw("not yet")
    elseif typeof(bc2)<:MatchedBC{2}
        Σ0R = kron(Iy,sparse([Nx],[Nx],[+αx_half⁻¹[Nx+1]/dx^2],Nx,Nx))
        Σ1R = kron(Iy,sparse([Nx],[Nx],[+1/dx],Nx,Nx))
        if isempty(bc2.in)
            fR = FΣ2D(dx,+1)
        else
            fR = FΣ2D(dx,-1)
        end
    else
        Σ0R = spzeros(N,N)
        Σ1R = spzeros(N,N)
        fR = FΣ2D(dx,0)
    end

    bc3 = domain.boundary.bcs[3]
    if typeof(bc3)<:DirichletBC{3}
        Σ0B = kron(sparse([1],[1],[-αy_half⁻¹[1]/dy^2],Ny,Ny),Ix)
        Σ1B = kron(sparse([1],[1],[+1/dy],Ny,Ny),Ix)
        fB = FΣ2D(dy,0)
    elseif typeof(bc3)<:NeumannBC{3}
        Σ0B = kron(sparse([1],[1],[+αy_half⁻¹[1]/dy^2],Ny,Ny),Ix)
        Σ1B = kron(sparse([1],[1],[-1/dy],Ny,Ny),Ix)
        fB = FΣ2D(dy,0)
    elseif typeof(bc3)<:FloquetBC{3}
        # Σ0L = sparse([li],[ri],[α_half⁻¹[ri+1]/dx^2],N,N)
        # Σ1L = sparse([li],[ri],[-1/dx],N,N)
        # fL = (ω,ka,_...) -> ffl(ka,a)
        throw("not yet")
    elseif typeof(bc3)<:MatchedBC{3}
        Σ0B = kron(sparse([1],[1],[+αy_half⁻¹[1]/dy^2],Ny,Ny),Ix)
        Σ1B = kron(sparse([1],[1],[-1/dy],Ny,Ny),Ix)
        if isempty(bc3.in)
            fB = FΣ2D(dy,+1)
        else
            fB = FΣ2D(dy,-1)
        end
    else
        Σ0B = spzeros(N,N)
        Σ1B = spzeros(N,N)
        fB = FΣ2D(dy,0)
    end

    bc4 = domain.boundary.bcs[4]
    if typeof(bc4)<:DirichletBC{4}
        Σ0T = kron(sparse([Ny],[Ny],[-αy_half⁻¹[Ny+1]/dy^2],Ny,Ny),Ix)
        Σ1T = kron(sparse([Ny],[Ny],[-1/dy],Ny,Ny),Ix)
        fT = FΣ2D(dy,0)
    elseif typeof(bc4)<:NeumannBC{4}
        Σ0T = kron(sparse([Ny],[Ny],[+αy_half⁻¹[Ny+1]/dy^2],Ny,Ny),Ix)
        Σ1T = kron(sparse([Ny],[Ny],[+1/dy],Ny,Ny),Ix)
        fT = FΣ2D(dy,0)
    elseif typeof(bc4)<:FloquetBC{4}
        # Σ0R = sparse([ri],[1],[-α_half⁻¹[li]/dx^2],N,N)
        # Σ1R = sparse([ri],[1],[+1/dx],N,N)
        # fR = (ω,ka,_...) -> ffr(ka,a)
        throw("not yet")
    elseif typeof(bc4)<:MatchedBC{4}
        Σ0T = kron(sparse([Ny],[Ny],[+αy_half⁻¹[Ny+1]/dy^2],Ny,Ny),Ix)
        Σ1T = kron(sparse([Ny],[Ny],[+1/dy],Ny,Ny),Ix)
        if isempty(bc4.in)
            fT = FΣ2D(dy,+1)
        else
            fT = FΣ2D(dy,-1)
        end
    else
        Σ0T = spzeros(N,N)
        Σ1T = spzeros(N,N)
        fT = FΣ2D(dy,0)
    end

    Σ2 = spzeros(N,N)

    fs = (fL,fR,fB,fT)

    return SelfEnergy{2,Symmetric,4,typeof(fs)}((-Σ0L,-Σ0R,-Σ0B,-Σ0T),(Σ1L,Σ1R,Σ1B,Σ1T),Σ2,fs)
end

struct FΣ2D <: Function
    dx::Float64
    sign::Int
end
function (fm::FΣ2D)(ω::Number,ka,ky=0,kz=0)
    if iszero(fm.sign)
        return complex(1.0,0)
    else
        k̂ = sqrt(ω^2 - ky^2 - kz^2)*fm.dx/2
        abs(k̂)<1 || throw(DomainError(fm.dx,"lattice spacing $(fm.dx) too large. for |ω|=$(abs(ω)), ky=$ky, kz=$kz, should be < $(abs(1/sqrt(ω^2 - ky^2 - kz^2)))"))
        kx = 2asin(k̂)/fm.dx
        return exp(fm.sign*1im*kx*fm.dx)
    end
end
function (fm::FΣ2D)(ω::Matrix,ka,ky,kz)
    if iszero(fm.sign)
        return one(convert(Matrix{ComplexF64},ω))
    else
        kx::Matrix{ComplexF64} = 2asin(sqrt(ω^2 - one(ω)*ky^2 - one(ω)*kz^2)*fm.dx/2)/fm.dx
        return exp(fm.sign*1im*kx*fm.dx)
    end
end

end

using .SelfEnergy2DSymmetric

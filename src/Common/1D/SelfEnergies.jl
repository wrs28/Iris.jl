function SelfEnergy(domains::NTuple{L,Domain{1}},surface,domain_index,α_half,a,nnm,nnp,indices,interior,domain_wall) where L
    li = findfirst(surface)
    ri = findlast(surface)
    # surface[1] || throw("first site should be an endpoint, something is amiss")
    # surface[end] || throw("last site should be an endpoint, something is amiss")
    sum(surface)>2 && throw("no point except first and last should be an endpoint, something is amiss")

    N = length(surface)
    α_half⁻¹ = 1 ./α_half

    dx1 = domains[domain_index[li]].lattice.dx
    bc1 = domains[domain_index[li]].boundary.bcs[1]

    if typeof(bc1)<:DirichletBC{1}
        Σ0L = sparse([li],[li],[-α_half⁻¹[li]/dx1^2],N,N)
        Σ1L = sparse([li],[li],[+1/dx1],N,N)
        fL = FM(dx1,0)
    elseif typeof(bc1)<:NeumannBC{1}
        Σ0L = sparse([li],[li],[+α_half⁻¹[li]/dx1^2],N,N)
        Σ1L = sparse([li],[li],[-1/dx1],N,N)
        fL = FM(dx1,0)
    elseif typeof(bc1)<:FloquetBC{1}
        # Σ0L = sparse([li],[ri],[α_half⁻¹[ri+1]/dx1^2],N,N)
        # Σ1L = sparse([li],[ri],[-1/dx1],N,N)
        # fL = (ω,ka,_...) -> ffl(ka,a)
        throw("not yet")
    elseif typeof(bc1)<:MatchedBC{1}
        Σ0L = sparse([li],[li],[+α_half⁻¹[li]/dx1^2],N,N)
        Σ1L = sparse([li],[li],[-1/dx1],N,N)
        if isempty(bc1.in)
            fL = FM(dx1,+1)
        else
            fL = FM(dx1,-1)
        end
    else
        Σ0L = spzeros(N,N)
        Σ1L = spzeros(N,N)
        fL = FM(dx1,0)
    end

    dx2 = domains[domain_index[ri]].lattice.dx
    bc2 = domains[domain_index[ri]].boundary.bcs[2]
    if typeof(bc2)<:DirichletBC{2}
        Σ0R = sparse([ri],[ri],[-α_half⁻¹[ri+1]/dx2^2],N,N)
        Σ1R = sparse([ri],[ri],[-1/dx2],N,N)
        fR = FM(dx2,0)
    elseif typeof(bc2)<:NeumannBC{2}
        Σ0R = sparse([ri],[ri],[+α_half⁻¹[ri+1]/dx2^2],N,N)
        Σ1R = sparse([ri],[ri],[+1/dx2],N,N)
        fR = FM(dx2,0)
    elseif typeof(bc2)<:FloquetBC{2}
        # Σ0R = sparse([ri],[1],[-α_half⁻¹[li]/dx2^2],N,N)
        # Σ1R = sparse([ri],[1],[+1/dx2],N,N)
        # fR = (ω,ka,_...) -> ffr(ka,a)
        throw("not yet")
    elseif typeof(bc2)<:MatchedBC{2}
        Σ0R = sparse([ri],[ri],[+α_half⁻¹[ri+1]/dx2^2],N,N)
        Σ1R = sparse([ri],[ri],[+1/dx2],N,N)
        if isempty(bc2.in)
            fR = FM(dx2,+1)
        else
            fR = FM(dx2,-1)
        end
    else
        Σ0R = spzeros(N,N)
        Σ1R = spzeros(N,N)
        fR = FM(dx1,0)
    end
    Σ2 = spzeros(N,N)

    fs = (fL,fR)

    return SelfEnergy{1,2,typeof(fs)}((-Σ0L,-Σ0R),(Σ1L,Σ1R),Σ2,(fL,fR))
end

# fdl(ω,_...) = one(ω)
# fnl(ω,_...) = one(ω)
# ffl(ka,a) = exp(-1im*ka*a)

# fdr(ω,_...) = one(ω)
# fnr(ω,_...) = one(ω)
# ffr(ka,a) = exp(+1im*ka*a)

# fmin(ω,ka,ky,kz,dx)  = exp(-1im*sqrt(ω^2-one(ω)*ky^2-one(ω)*kz^2)*dx)
# fmout(ω,ka,ky,kz,dx) = exp(+1im*sqrt(ω^2-one(ω)*ky^2-one(ω)*kz^2)*dx)


struct FM
    dx::Float64
    sign::Int
end
function (fm::FM)(ω,ka,ky,kz)
    if iszero(fm.sign)
        complex(1,0)
    else
        exp(fm.sign*1im*sqrt(ω^2-one(ω)*ky^2-one(ω)*kz^2)*fm.dx)
    end
end

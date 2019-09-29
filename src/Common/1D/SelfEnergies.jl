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
        fL = fdl
    elseif typeof(bc1)<:NeumannBC{1}
        Σ0L = sparse([li],[li],[+α_half⁻¹[li]/dx1^2],N,N)
        Σ1L = sparse([li],[li],[-1/dx1],N,N)
        fL = fnl
    elseif typeof(bc1)<:FloquetBC{1}
        Σ0L = sparse([li],[ri],[α_half⁻¹[ri+1]/dx1^2],N,N)
        Σ1L = sparse([li],[ri],[-1/dx1],N,N)
        fL = (ω,ka,_...) -> ffl(ka,a)
    elseif typeof(bc1)<:MatchedBC{1}
        Σ0L = sparse([li],[li],[+α_half⁻¹[li]/dx1^2],N,N)
        Σ1L = sparse([li],[li],[-1/dx1],N,N)
        if isempty(bc1.in)
            fL = (ω...) -> fmout(ω...,dx1)
        else
            fL = (ω...) -> fmin(ω...,dx1)
        end
    else
        Σ0L = spzeros(N,N)
        Σ1L = spzeros(N,N)
        fL = fnl
    end

    dx2 = domains[domain_index[ri]].lattice.dx
    bc2 = domains[domain_index[ri]].boundary.bcs[2]
    if typeof(bc2)<:DirichletBC{2}
        Σ0R = sparse([ri],[ri],[-α_half⁻¹[ri+1]/dx2^2],N,N)
        Σ1R = sparse([ri],[ri],[-1/dx2],N,N)
        fR = fdr
    elseif typeof(bc2)<:NeumannBC{2}
        Σ0R = sparse([ri],[ri],[+α_half⁻¹[ri+1]/dx2^2],N,N)
        Σ1R = sparse([ri],[ri],[+1/dx2],N,N)
        fR = fnr
    elseif typeof(bc2)<:FloquetBC{2}
        Σ0R = sparse([ri],[1],[-α_half⁻¹[li]/dx2^2],N,N)
        Σ1R = sparse([ri],[1],[+1/dx2],N,N)
        fR = (ω,ka,_...) -> ffr(ka,a)
    elseif typeof(bc2)<:MatchedBC{2}
        Σ0R = sparse([ri],[ri],[+α_half⁻¹[ri+1]/dx2^2],N,N)
        Σ1R = sparse([ri],[ri],[+1/dx2],N,N)
        if isempty(bc2.in)
            fR = (ω...) -> fmout(ω...,dx2)
        else
            fR = (ω...) -> fmin(ω...,dx2)
        end
    else
        Σ0R = spzeros(N,N)
        Σ1R = spzeros(N,N)
    end
    Σ2 = spzeros(N,N)

    fs = (fL,fR)

    return SelfEnergy{1,2,typeof(fs)}((Σ0L,Σ0R),(Σ1L,Σ1R),Σ2,(fL,fR))
end

fdl(ω,_...) = one(ω)
fnl(ω,_...) = one(ω)
ffl(ka,a) = exp(-1im*ka*a)

fdr(ω,_...) = one(ω)
fnr(ω,_...) = one(ω)
ffr(ka,a) = exp(+1im*ka*a)

fmin(ω,ka,ky,kz,dx)  = exp(-1im*sqrt(ω^2-one(ω)*ky^2-one(ω)*kz^2)*dx)
fmout(ω,ka,ky,kz,dx) = exp(+1im*sqrt(ω^2-one(ω)*ky^2-one(ω)*kz^2)*dx)

module Scattering1D

using LinearAlgebra
using SparseArrays
using ..Common
import ..Common.EQUIVALENT_SOURCE_RELATIVE_CUTOFF
import ..Common.AbstractComplexBL
import ..LinearScatteringProblem
import ..NonLinearScatteringProblem
import ..HelmholtzLS
import ..HelmholtzNLS
import ..MaxwellLS
import ..MaxwellNLS
import ..EquivalentSource
import ..ScatteringSolution
import ..ψ_to_x


Common.update!(mls::LinearScatteringProblem{1}, aL::Number, aR::Number, args...) = update!(mls, [aL,aR], args...)

function Common.update!(mls::LinearScatteringProblem{1}, a::Vector, args...)
    if !isempty(args)
        mls.ω = args[1]
        mls.operator(mls.ω)
    end
    copyto!(mls.equivalent_source.a, a)
    fill!(mls.equivalent_source.field,0)
    mls.j[:] = mls.equivalent_source(mls.operator.D², mls.operator.ky, mls.operator.kz)
    return nothing
end


"""
    HelmholtzLS(::Simulation{1}, ω, aL, [aR=0 ;  start=-Inf, stop=Inf, ky=0, kz=0]) -> ls
"""
HelmholtzLS(sim::Simulation{1}, ω::Real, aL::Number, aR::Number=0; kwargs...) = HelmholtzLS(sim,ω,[aL,aR];kwargs...)

"""
    HelmholtzLS(::Simulation{1}, ω, [amplitudes::Vector ;  start=-Inf, stop=Inf, ky=0, kz=0]) -> ls
"""
function HelmholtzLS(sim::Simulation{1}, ω::Real, a::Vector; start::Real=-Inf, stop::Real=Inf)
    M = Helmholtz(sim)
    M(ω) # evaluate at ω
    es = EquivalentSource{1}(sim, ω, a; start=start, stop=stop)
    j = es(M.D²,0,0)
    sol = ScatteringSolution{1}(sim,ω,a)
    return LinearScatteringProblem{1,1,typeof(M),typeof(es)}(M,es,j,Ref(false),Ref(false),sol)
end


"""
    MaxwellLS(::Simulation{1}, ω, aL, [aR=0 ;  start=-Inf, stop=Inf, ky=0, kz=0]) -> ls
"""
MaxwellLS(sim::Simulation{1}, ω::Real, aL::Number, aR::Number=0; kwargs...) = MaxwellLS(sim,ω,[aL,aR];kwargs...)

"""
    MaxwellLS(::Simulation{1}, ω, [amplitudes::Vector ;  start=-Inf, stop=Inf, ky=0, kz=0]) -> ls
"""
function MaxwellLS(sim::Simulation{1}, ω::Real, a::Vector; start::Real=-Inf, stop::Real=Inf, ky::Number=0, kz::Number=0)
    M = Maxwell(sim; ky=ky,kz=kz)
    M(ω) # evaluate at ω
    es = EquivalentSource{3}(sim, ω, a; start=start, stop=stop, ky=ky, kz=kz)
    j = es(M.D²,ky,kz)
    sol = ScatteringSolution{3}(sim,ω,a)
    return LinearScatteringProblem{1,3,typeof(M),typeof(es)}(M,es,j,Ref(false),Ref(false),sol)
end


"""
    HelmholtzNLS(::Simulation{1}, ω, aL, [aR=0, lupack=$DEFAULT_LUPACK;  start=-Inf, stop=Inf, ky=0, kz=0]) -> ls
"""
HelmholtzNLS(sim::Simulation{1}, ω::Real, aL::Number, aR::Number=0, lupack::AbstractLUPACK=DEFAULT_LUPACK; kwargs...) = HelmholtzNLS(sim,ω,[aL,aR]; kwargs...)

"""
    HelmholtzNLS(::Simulation{1}, ω, [amplitudes::Vector, lupack=$DEFAULT_LUPACK;  start=-Inf, stop=Inf, ky=0, kz=0]) -> ls
"""
function HelmholtzNLS(sim::Simulation{1}, ω::Real, a::Vector, lupack::AbstractLUPACK=DEFAULT_LUPACK; start::Real=-Inf, stop::Real=Inf)
    ls = HelmholtzLS(sim, ω, a; start=start, stop=stop)
    ψ = ScalarField(sim,1)
    res = similar(ψ)
    x = ψ_to_x(ψ)
    sol = ScatteringSolution{1}(sim,ω,a)
    return NonLinearScatteringProblem{1,1,typeof(ls),typeof(lupack)}(ls,ψ,res,x,Ref(false),lupack)
end

"""
    MaxwellNLS(::Simulation{1}, ω, aL, [aR=0, lupack=$DEFAULT_LUPACK;  start=-Inf, stop=Inf, ky=0, kz=0]) -> nls
"""
MaxwellNLS(sim::Simulation{1},ω::Real,aL::Number,aR::Number=0, lupack::AbstractLUPACK=DEFAULT_LUPACK; kwargs...) = MaxwellNLS(sim,ω,[aL,aR]; kwargs...)

"""
    MaxwellNLS(::Simulation{1}, ω, [amplitudes::Vector, lupack=$DEFAULT_LUPACK;  start=-Inf, stop=Inf, ky=0, kz=0]) -> nls
"""
function MaxwellNLS(sim::Simulation{1}, ω::Real, a::Vector, lupack::AbstractLUPACK=DEFAULT_LUPACK; start::Real=-Inf, stop::Real=Inf, ky::Number=0, kz::Number=0)
    ls = MaxwellLS(sim, ω, a; start=start, stop=stop, ky=ky, kz=kz)
    ψ = ElectricField(sim,1)
    res = similar(ψ)
    x = ψ_to_x(ψ)
    sol = ScatteringSolution{3}(sim,ω,a)
    return NonLinearScatteringProblem{1,3,typeof(ls),typeof(lupack)}(ls,ψ,res,x,Ref(false),lupack)
end

################################################################################
# EQUIVALENT SOURCE BLOCK FOR DIM=1

# constructor
function EquivalentSource{M}(sim::Simulation{1}, ω::Number, a::Vector; start::Number=-Inf,stop::Number=Inf, ky::Number=0, kz::Number=0) where M
    field = VectorField{M}(sim,1)
    bc1, bc2 = sim.boundary.bcs
    bl1, bl2 = sim.boundary.bls
    if typeof(bl1)<:noBL{1}
        start_in = start
        start_out = start
    else
        start_in  = max(start, sim.boundary.bls[1].start + 5sim.dx + 3sim.dx)
        start_out = max(start, sim.boundary.bls[1].start + 5sim.dx - 3sim.dx)
    end
    if typeof(bl2)<:noBL{2}
        stop_in = stop
        stop_out = stop
    else
        stop_in  = min(stop, sim.boundary.bls[2].start - 5sim.dx - 3sim.dx)
        stop_out = min(stop, sim.boundary.bls[2].start - 5sim.dx + 3sim.dx)
    end
    incoming_mask = Interval(max(field.start,start_in),min(field.stop,stop_in))
    outgoing_mask = Interval(max(field.start,start_out),min(field.stop,stop_out))
    kx = 2asin(sqrt(ω^2 - ky^2 - kz^2)*sim.dx/2)/sim.dx
    channelflux = fill(hypot(kx,ky,kz),2)
    return EquivalentSource(incoming_mask,outgoing_mask,field,sim,ω,a,channelflux)
end

# generate source
function (es::EquivalentSource{1,MCOMPONENTS})(ky::Number=0,kz::Number=0) where MCOMPONENTS
    if MCOMPONENTS==1
        M = Helmholtz(es.sim)
    elseif MCOMPONENTS==3
        M = Maxwell(es.sim;ky=ky,kz=kz)
    end
    M(es.ω)
    return es(M.D²,ky,kz)
end

# generate source, given differential operator
@inline function (es::EquivalentSource{1,M,TW,TA})(D²,ky::Number,kz::Number) where {M,TA<:AbstractVecOrMat,TW<:Number}
    n = length(es.sim)
    ω, a = es.ω, es.a
    ω² = ω^2
    kx = 2asin(sqrt(ω² - ky^2 - kz^2)*es.sim.dx/2)/es.sim.dx

    bl1, bl2 = es.sim.boundary.bls
    bc1, bc2 = es.sim.boundary.bcs

    # left  = es.sim.x[1] - es.sim.dx/2
    # right = es.sim.x[end] + es.sim.dx/2
    if typeof(bl1)<:AbstractComplexBL{1} || typeof(bc1)<:MatchedBC{1}
        left  = es.field.start
    else
        left  = es.sim.x[1] - es.sim.dx/2
    end
    if typeof(bl2)<:AbstractComplexBL{2} || typeof(bc2)<:MatchedBC{2}
        right = es.field.stop
    else
        right = es.sim.x[end] + es.sim.dx/2
    end

    temp = Vector{ComplexF64}(undef,M*n)
    rows = Vector{Int}(undef,M*n)
    count = 0
    # @fastmath @inbounds
    for i ∈ eachindex(es.field)
        j = mod1(i,n)
        if (i-1)÷n < 2
            if M==3
                β = (i-1)÷n==0 ? (ky/ω) : 1.0
            elseif M==1
                β = 1
            end
            incoming = es.incoming_mask(es.sim.x[j])
            outgoing = es.outgoing_mask(es.sim.x[j])
            es.field[i] = 0
            if incoming
                if typeof(bl2)<:AbstractComplexBL{2} || typeof(bc2)<:MatchedBC{2}
                    nothing
                elseif typeof(bc2)<:Union{DirichletBC{2},NeumannBC{2}}
                    es.field[i] += a[1]*β*cis(kx*(es.field.pos[j].x-right.x))/sqrt(es.channelflux[1])
                    iszero(a[2]) || throw("cannot be incident from $(typeof(bc2))")
                elseif typeof(bc2)<:FloquetBC{2}
                    throw("Floquet{2} incompatible with scattering")
                end
                if typeof(bl1)<:AbstractComplexBL{1} || typeof(bc1)<:MatchedBC{1}
                    nothing
                elseif typeof(bc1)<:Union{DirichletBC{1},NeumannBC{1}}
                    es.field[i] += a[2]*β*cis(-kx*(es.field.pos[j].x-left.x))/sqrt(es.channelflux[2])
                    iszero(a[1]) || throw("cannot be incident from $(typeof(bc1))")
                elseif typeof(bc1)<:FloquetBC{1}
                    throw("Floquet{1} incompatible with scattering")
                end
            end
            if outgoing
                if typeof(bl2)<:AbstractComplexBL{2} || typeof(bc2)<:MatchedBC{2}
                    es.field[i] = a[1]*β*cis(kx*(es.field.pos[j].x-left.x))/sqrt(es.channelflux[1])
                elseif typeof(bc2)<:DirichletBC{2}
                    es.field[i] += -a[1]*β*cis(-kx*(es.field.pos[j].x-right.x))/sqrt(es.channelflux[1])
                elseif typeof(bc2)<:NeumannBC{2}
                    es.field[i] += a[1]*β*cos(-kx*(es.field.pos[j].x-right.x))/sqrt(es.channelflux[1])
                end
                if typeof(bl1)<:AbstractComplexBL{1} || typeof(bc1)<:MatchedBC{1}
                    es.field[i] += a[2]*β*cis(-kx*(es.field.pos[j].x-right.x))/sqrt(es.channelflux[2])
                elseif typeof(bc1)<:DirichletBC{1}
                    es.field[i] += -a[2]*β*cis(kx*(es.field.pos[j].x-left.x))/sqrt(es.channelflux[2])
                elseif typeof(bc1)<:NeumannBC{1}
                    es.field[i] += a[2]*β*cis(kx*(es.field.pos[j].x-left.x))/sqrt(es.channelflux[2])
                end
            end
            if incoming || outgoing
                count += 1
                temp[count] = es.field[i]
                rows[count] = i
            end
        end
    end

    if M==3
        for i ∈ diagind(D²) D²[i] -= ω² end
    elseif M==1
        for i ∈ diagind(D²) D²[i] += ω² end
    end
    j = D²*sparsevec(rows[1:count],temp[1:count],M*n)
    if M==3
        for i ∈ diagind(D²) D²[i] += ω² end
    elseif M==1
        for i ∈ diagind(D²) D²[i] -= ω² end
    end

    mx = maximum(abs,j)
    return droptol!(j,mx*EQUIVALENT_SOURCE_RELATIVE_CUTOFF)
end

end # module


using .Scattering1D

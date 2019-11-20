"""
    scattering(::Simulation{1}, ω, left_amplitude, [aright_amplitude=0; lupack=$(DEFAULT_LUPACK), start=-Inf, stop=Inf, ky=0, kz=0]) -> sol
"""
scattering(sim::Simulation{1}, ω::Real, aL::Number, aR::Number=0; kwargs...) =
    scattering(sim,ω,[aL,aR]; kwargs...)

"""
    scattering(::Simulation{1}, ω, amplitudes::Vector;[lupack=$(DEFAULT_LUPACK), start=-Inf, stop=Inf, ky=0, kz=0]) -> sol
"""
function scattering(
            sim::Simulation{1},
            ω::Real,
            a::AbstractVecOrMat;
            lupack::AbstractLUPACK=DEFAULT_LUPACK,
            kwargs...)

    mls = Maxwell_LS(sim, ω, a; kwargs...)
    scattering!(mls, lupack)
    return mls.solution
end

"""
    maxwell_ls(::Simulation{1}, ω, aL, [aR=0 ;  start=-Inf, stop=Inf, ky=0, kz=0]) -> nls
    maxwell_ls(::Simulation{1}, ω, [amplitudes::Vector ;  start=-Inf, stop=Inf, ky=0, kz=0]) -> nls
"""
maxwell_ls
Maxwell_LS(sim::Simulation{1}, ω::Real, aL::Number, aR::Number=0; kwargs...) =
    Maxwell_LS(sim,ω,[aL,aR];kwargs...)
function Maxwell_LS(sim::Simulation{1},ω::Real,a::Vector;
            start::Real=-Inf, stop::Real=Inf, ky::Number=0, kz::Number=0)
    M = maxwell(sim; ky=ky,kz=kz)
    M(ω) # evaluate at ω
    es = EquivalentSource(sim, ω, a; start=start, stop=stop, ky=ky, kz=kz)
    j = es(M.D²,ky,kz)
    sol = ScatteringSolution(sim,ω,a)
    return Maxwell_LS{1,typeof(M),typeof(es)}(M,es,j,Ref(false),Ref(false),sol)
end


Common.update!(mls::Maxwell_LS{1}, aL::Number,aR::Number) = update!(mls,[aL,aR])
function Common.update!(mls::Maxwell_LS{1}, a::Vector)
    mls.es.a[:] = a
    mls.j[:] = es(mls.M.D², mls.maxwell.ky, mls.maxwell.kz)
    return nothing
end


################################################################################
# nonlinear scattering convenience wrappers

scattering_nl(sim::Simulation{1},ω::Real,aL::Number,aR::Number=0; kwargs...) =
    scattering_nl(sim,ω,[aL,aR];kwargs...)

function scattering_nl(
    sim::Simulation{1},
    ω::Real,
    a::AbstractVecOrMat;
    start::Real = -Inf,
    stop::Real = Inf,
    ky::Number=0,
    kz::Number=0,
    init::ElectricField{1}=scattering(sim,ω,a;start=start,stop=stop,ky=ky,kz=kz).tot,
    kwargs...)

    sol = ScatteringSolution(sim,ω,a)
    scattering_nl!(sol, sim; start=start, stop=stop, ky=ky, kz=kz, init=init, kwargs...)
    return sol
end


function scattering_nl!(
            sol::ScatteringSolution{1},
            sim::Simulation{1};
            start::Real=-Inf,
            stop::Real=Inf,
            ky::Number=0,
            kz::Number=0,
            init::ElectricField{1},
            verbose::Bool=false,
            show_trace::Bool=verbose,
            lupack::AbstractLUPACK=DEFAULT_LUPACK,
            kwargs...)

    nls = maxwell_nls(sim, sol.ω, sol.a, lupack; start=start, stop=stop, ky=ky, kz=kz)
    results = nlsolve(nls, init; show_trace=show_trace, kwargs...)
    x_to_ψ!(nls,results.zero)
    (results.f_converged || results.x_converged) || printstyled("beware: no convergence at ω=$(nls.ω), a=$(nls.a)\n", color=PRINTED_COLOR_BAD)
    field = nls.equivalent_source.field
    for i ∈ eachindex(nls.equivalent_source.field)
        k = mod1(i,length(sim))
        sol.total[i] = nls.ψ[i]
        sol.incident[i] = field[i]*nls.equivalent_source.mask(field.pos[k])
        sol.scattered[i] = sol.total[i] - sol.incident[i]
    end
    return nothing
end

# constructors for DIM=1
"""
    maxwell_nls(::Simulation{1}, ω, aL, [aR=0, lupack=$DEFAULT_LUPACK;  start=-Inf, stop=Inf, ky=0, kz=0]) -> nls
    maxwell_nls(::Simulation{1}, ω, [amplitudes::Vector, lupack=$DEFAULT_LUPACK;  start=-Inf, stop=Inf, ky=0, kz=0]) -> nls
"""
maxwell_nls
Maxwell_NLS(sim::Simulation{1},ω::Real,aL::Number,aR::Number=0, lupack::AbstractLUPACK=DEFAULT_LUPACK; kwargs...) =
    Maxwell_NLS(sim,ω,[aL,aR]; kwargs...)
function Maxwell_NLS(sim::Simulation{1},ω::Real,a::Vector,lupack::AbstractLUPACK=DEFAULT_LUPACK; start::Real=-Inf, stop::Real=Inf, ky::Number=0, kz::Number=0)
    ls = Maxwell_LS(sim, ω, a; start=start, stop=stop, ky=ky, kz=kz)
    ψ = ElectricField(sim,1)
    res = similar(ψ)
    x = ψ_to_x(ψ)
    sol = ScatteringSolution(sim,ω,a)
    return Maxwell_NLS{1,typeof(ls),typeof(lupack)}(ls,ψ,res,x,Ref(false),lupack)
end

################################################################################
# EQUIVALENT SOURCE BLOCK FOR DIM=1

# constructor
function EquivalentSource(sim::Simulation{1}, ω::Number, a ;start::Number=-Inf,stop::Number=Inf,ky::Number=0,kz::Number=0)
    mask = Interval(start,stop)
    field = ElectricField(sim,zeros(ComplexF64,3length(sim)))
    return EquivalentSource(mask,field,sim,ω,a)
end

# generate source
function (es::EquivalentSource{1})(ky=0,kz=0)
    M = maxwell(es.sim;ky=ky,kz=kz)
    M(es.ω)
    return es(M.D²,ky,kz)
end

# generate source, given differential operator
@inline function (es::EquivalentSource{1,TW,TA})(D²,ky,kz) where {TA<:AbstractVecOrMat,TW<:Number}
    n = length(es.sim)
    ω, a = es.ω, es.a
    ω² = ω^2
    kx = 2asin(sqrt(ω² - ky^2 - kz^2)*es.sim.dx/2)/es.sim.dx

    left = findmin(map(s->s.x,es.sim.x))[1]-es.sim.dx/2
    right = findmax(map(s->s.x,es.sim.x))[1]+es.sim.dx/2

    temp = Vector{ComplexF64}(undef,3n)
    rows = Vector{Int}(undef,3n)
    count = 0
    @fastmath @inbounds for i ∈ eachindex(es.field)
        j = mod1(i,n)
        if (i-1)÷n == 0 && es.mask(es.sim.x[j])
            count += 1
            es.field[i] = a[1]*(ky/ω)*cis(kx*(es.field.pos[j].x-left))/sqrt(ω)
            es.field[i] += a[2]*(ky/ω)*cis(-kx*(es.field.pos[j].x-right))/sqrt(ω)
            temp[count] = es.field[i]
            rows[count] = i
        elseif (i-1)÷n == 1 && es.mask(es.field.pos[j])
            count += 1
            es.field[i] = a[1]*cis(kx*(es.field.pos[j].x-left))/sqrt(ω)
            es.field[i] += a[2]*cis(-kx*(es.field.pos[j].x-right))/sqrt(ω)
            temp[count] = es.field[i]
            rows[count] = i
        end
    end

    for i ∈ diagind(D²) D²[i] -= ω² end
    j = D²*sparsevec(rows[1:count],temp[1:count],3n)
    for i ∈ diagind(D²) D²[i] += ω² end

    mx = maximum(abs,j)
    return droptol!(j,mx*EQUIVALENT_SOURCE_RELATIVE_CUTOFF)
end


################################################################################

include("Plotting1D.jl")

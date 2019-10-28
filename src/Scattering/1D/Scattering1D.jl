function scattering(
            sim::Simulation{1},
            ω::Real,
            lupack::Type{T}=DEFAULT_LUPACK;
            start::Real=-Inf,
            stop::Real=Inf,
            ky::Number=0,
            kz::Number=0) where T<:AbstractLUPACK

    n = 3length(sim)
    sol = ScatteringSolution(sim.x, Vector{Float64}(undef,n),Vector{Float64}(undef,n),Vector{Float64}(undef,n), ω)
    scattering!(sol, sim, ω, lupack; start=start, stop=stop, ky=ky, kz=kz)
    return sol
end

function scattering!(
            sol::ScatteringSolution{1},
            sim::Simulation{1},
            ω::Real,
            lupack::Type{T};
            start::Real=-Inf,
            stop::Real=Inf,
            ky::Number=0,
            kz::Number=0) where T<:AbstractLUPACK

    n = 3length(sim)
    sc = EquivalentSource(sim; start=start, stop=stop, ky=ky, kz=kz)
    A,B = maxwell_lep(sim;ky=ky,kz=kz)(ω)
    j = sc(ω,A,ky,kz)
    alu = lu(A-lmul!(ω^2,B),lupack)
    if lupack<:PSolver
        ldiv!(sol.total.val, alu, Matrix(reshape(j,n,1)))
    else
        ldiv!(sol.total.val, alu, j)
    end
    for i ∈ eachindex(sc.field)
        sol.incident[i] = sc.field[i]
        sol.scattered[i] = sol.total[i] - sol.incident[i]
    end
    return nothing
end


function scattering_nl(
    sim::Simulation{1},
    ω::Real;
    start::Real=-Inf,
    stop::Real=Inf,
    ky::Number=0,
    kz::Number=0,
    ψ_init::ElectricField{1}=scattering(sim,ω;start=start,stop=stop,ky=ky,kz=kz).tot,
    verbose::Bool=false,
    show_trace::Bool=verbose,
    kwargs...
    )

    n = 3length(sim)
    sol = ScatteringSolution(sim.x,Array{Float64}(undef,n),Array{Float64}(undef,n),Array{Float64}(undef,n),ω)
    scattering_nl!(sol,sim,ω;start=start,stop=stop,ky=ky,kz=kz,ψ_init=ψ_init,show_trace=show_trace,kwargs...)
    return sol
end


function scattering_nl!(
            sol::ScatteringSolution{1},
            sim::Simulation{1},
            ω::Real;
            start::Real=-Inf,
            stop::Real=Inf,
            ky::Number=0,
            kz::Number=0,
            ψ_init::ElectricField{1},
            verbose::Bool=false,
            show_trace::Bool=verbose,
            kwargs...)

    nls = maxwell_nls(sim,ω;start=start,stop=stop,ky=ky,kz=kz)
    results = nlsolve(nls, ψ_init; show_trace=show_trace,kwargs...)
    x_to_ψ!(nls,results.zero)
    (results.f_converged || results.x_converged) || printstyled("beware: no convergence",color=PRINTED_COLOR_BAD)
    field = nls.equivalent_source.field
    for i ∈ eachindex(nls.equivalent_source.field)
        sol.total[i] = nls.ψ[i]
        sol.incident[i] = field[i]
        sol.scattered[i] = sol.total[i] - sol.incident[i]
    end
    return nothing
end


function Maxwell_NLS(sim::Simulation{1},ω::Real; start::Real=-Inf, stop::Real=Inf, ky::Number=0, kz::Number=0)
    M = maxwell(sim;ky=ky,kz=kz)
    es = EquivalentSource(sim;start=start,stop=stop,ky=ky,kz=kz)
    j = es(ω)
    ψ = ElectricField(sim.x,1)
    res = similar(ψ)
    x = ψ_to_x(ψ)
    return Maxwell_NLS{1,typeof(M),typeof(es)}(ω,M,es,j,ψ,res,x,Ref(false),Ref(false))
end

################################################################################

function EquivalentSource(sim::Simulation{1};start::Number=-Inf,stop::Number=Inf,ky::Number=0,kz::Number=0)
    n = length(sim)
    mask = Interval(start,stop)
    field = ElectricField(sim.x,zeros(ComplexF64,3n))
    return EquivalentSource(mask,field,sim)
end

function (es::EquivalentSource{1})(ω::Number;ky::Number=0,kz::Number=0)
    n = length(es.sim)
    D² = maxwell_lep(es.sim;ky=ky,kz=kz)(ω)[1]
    return es(ω,D²,ky,kz)
end
function (es::EquivalentSource{1})(ω::Number,D²,ky,kz)
    n = length(es.sim)
    ω² = ω^2
    kx = 2asin(sqrt(ω² - ky^2 - kz^2)*es.sim.dx/2)/es.sim.dx
    for i ∈ diagind(D²) D²[i] -= ω² end
    temp = Vector{ComplexF64}(undef,3n)
    rows = Vector{Int}(undef,3n)
    count = 0
    for i ∈ eachindex(es.field)
        j = mod1(i,n)
        if (i-1)÷n == 0 && es.mask(es.sim.x[j])
            count += 1
            es.field[i] = (ky/ω)*exp(1im*kx*es.field.pos[j].x)/sqrt(ω)
            temp[count] = es.field[i]
            rows[count] = i
        elseif (i-1)÷n == 1 && es.mask(es.field.pos[j])
            count += 1
            es.field[i] = exp(1im*kx*es.field.pos[j].x)/sqrt(ω)
            temp[count] = es.field[i]
            rows[count] = i
        end
    end
    j = D²*sparsevec(rows[1:count],temp[1:count],3n)
    for i ∈ diagind(D²) D²[i] += ω² end

    mx = maximum(abs,j)
    for i ∈ eachindex(j.nzval)
        abs(j.nzval[i])/mx < EQUIVALENT_SOURCE_RELATIVE_CUTOFF ? j.nzval[i] = 0 : nothing
    end
    return dropzeros!(j)
end

include("Plotting1D.jl")

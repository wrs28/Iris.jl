# TODO: SALTbootstrap
# TODO: convenience Simulation modifiers to take all bc's to outgoing.
"""
    module Lasing

For solution of the nonlinear SALT equation.
Exports `HelmholtzSALT` and `MaxwellSALT` constructor and extends methods of
`NLsolve`, which must be separately imported.
"""
module Lasing

export HelmholtzSALT
export MaxwellSALT

using ..Common
using ..Spectral
using Formatting
using LinearAlgebra
using NLsolve
using RecipesBase
using SparseArrays

import ..Common.INDEX_OFFSET

################################################################################
# SALT object definitions and constructors

# works for Helmholtz and Maxwell, field `operator` is accessed through
# properties `helmholtz` or `maxwell`
struct SALTProblem{NMODES,M,TO,DIM}
    operator::TO
    ωs::Vector{Float64}
    ψs::VectorField{DIM,M}
    x::Vector{Float64}
    n::Int
    m::Int
    res::Vector{ComplexF64}
    index::Int
    solved::Ref{Bool}
    converged::Ref{Bool}
end

# main SALTProblem constructor
function SALTProblem{NMODES,MCOMPONENTS}(sim::Simulation{DIM}) where {DIM,NMODES,MCOMPONENTS}
    M = Helmholtz(sim; m=NMODES)
    n = length(sim)
    ωs = Vector{Float64}(undef, NMODES)
    ψs = VectorField{MCOMPONENTS}(sim.x, NMODES)
    n÷2>INDEX_OFFSET || throw("INDEX_OFFSET=$(INDEX_OFFSET÷MCOMPONENTS) requires number of sites in ψs_init ($n) to be greater than $(2INDEX_OFFSET÷MCOMPONENTS)")
    x = Vector{Float64}(undef,2n*NMODES)
    res = Vector{ComplexF64}(undef,n)
    index = n÷2+INDEX_OFFSET
    return SALTProblem{NMODES,MCOMPONENTS,typeof(M),DIM}(M,ωs,ψs,x,n,NMODES,res,index,Ref(false),Ref(false))
end


HelmholtzSALT{NMODES} = SALTProblem{NMODES,1}
"""
    HelmholtzSALT(sim::Simulation, m::Integer) -> salt

SALT object with `m` fields for Helmholtz problem
"""
HelmholtzSALT(sim::Simulation,m::Integer) = HelmholtzSALT{m}(sim)


MaxwellSALT{NMODES} = SALTProblem{NMODES,3}
"""
    MaxwellSALT(sim::Simulation, m::Integer) -> salt

SALT object with `m` fields for Maxwell problem
"""
MaxwellSALT(sim::Simulation,m::Integer) = MaxwellSALT{m}(sim)


function Base.getproperty(s::SALTProblem{N,M}, sym::Symbol) where {N,M}
    if Base.sym_in(sym,(:ω,:freq,:Freq,:frequencies,:frequency,:Frequencies,:Frequency,:omega,:omegas))
        return getfield(s,:ωs)
    elseif Base.sym_in(sym,(:ψ,:field,:Field,:fields,:Fields,:psi,:psis))
        return getfield(s,:ψs)
    elseif Base.sym_in(sym,(:sol,:Sol,:solution,:Solution))
        return (getfield(s,:ωs),getfield(s,:ψs))
    elseif Base.sym_in(sym,(:hemlholtz,:maxwell))
        return getfield(s,:operator)
    elseif sym == :simulation
        return getfield(getfield(s,:operator),:simulation)
    elseif Base.sym_in(sym,(:helmholtz,:maxwell))
        return getfield(s,:operator)
    else
        return getfield(s,sym)
    end
end

function Base.propertynames(::SALTProblem,private=false)
    if private
        return fieldnames(SALTProblem)
    else
        return (:ωs,:ψs,:solution,:simulation)
    end
end

################################################################################
# NLsolve extension

fnames = (:nlsolve,:fixedpoint)

for fn ∈ fnames

    @eval NLsolve.$(fn)(ms::SALTProblem; kwargs...) = $(fn)(ms,ms.ωs,ms.ψs; kwargs...)

    @eval NLsolve.$(fn)(ms::SALTProblem, ωs_init::Real, ψs_init::VectorField; kwargs...) = $(fn)(ms,[ωs_init],ψs_init; kwargs...)

    @eval begin function NLsolve.$(fn)(ms::SALTProblem{N,M}, ωs_init::Vector, ψs_init::VectorField; kwargs...) where {N,M}
            n,m = size(ψs_init)
            N==m || throw("number of modes in ψs_init ($m) must be the same as given in Maxwell_SALT ($N)")
            n==M*length(ms.simulation) || throw("number of sites in ψs_init ($(n÷M)) must be the same as given in simulation ($(length(ms.simulation)))")
            n÷2>INDEX_OFFSET || throw("INDEX_OFFSET=$(INDEX_OFFSET) requires number of sites in ψs_init ($(n÷M)) to be greater than $(2INDEX_OFFSET÷M)")
            ωψ_to_x!(ms.x,ωs_init,ψs_init)
            x_to_ωψ!(ms)
            F0 = similar(ms.x)
            J0 = initialize_J(ms)
            df = OnceDifferentiable(ms,ms,ms,ms.x,F0,J0)
            results = $(fn)(df,ms.x; kwargs...)
            # results = $(fn)(ms,ms.x; kwargs...)
            x_to_ωψ!(ms,results.zero)
            ms.solved[] = true
            ms.converged[] = results.f_converged || results.x_converged
            return results
        end
    end
end


"""
    nlsolve(salt; kwargs...) -> results

    nlsolve(salt, init_ωs, init_ψs; kwargs...) -> results

solve nonlinear SALT equation where `salt` is a `HelmholtzSALT` or `MaxwellSALT`.
Results also stored `salt`.

`init_ψs` must be a `ScalarField` (for Helmholtz) or `ElectricField` (for Maxwell)
"""
nlsolve

"""
    fixedpoint(salt; kwargs...) -> results

    fixedpoint(salt, init_ωs, init_ψs; kwargs...) -> results

Convenience wrapper for `NLsolve`'s `fixedpoint` wrapper.
Solves nonlinear SALT equation where `scpa` is a `HelmholtzSALT` or `MaxwellSALT`.
Results also stored `scpa`.

`init_ψs` must be a `ScalarField` (for Helmholtz) or `ElectricField` (for Maxwell)
"""
fixedpoint

################################################################################
# General nonlinear objective + jacobian
(s::SALTProblem)(F::Vector,x::Vector) = s(F,nothing,x)
(s::SALTProblem)(J::AbstractMatrix,x::Vector) = s(nothing,J,x)

################################################################################
# Helmholtz nonlinear objective + jacobian

# Helmholtz Jacobians/Objective
@inline function (ms::HelmholtzSALT)(F,J,x::Vector)
    x_to_ωψ!(ms,x)
    n,m = ms.n,ms.m
    nm = n*m
    index = ms.index

    if !isnothing(J) ψ = ms.ψs end

    ψμ = ms.ψs[:,1]
    @inbounds for μ ∈ 1:m
        Aμ = ms.helmholtz(ms.ωs[μ],ms.ωs,ms.ψs)
        @inbounds for i ∈ 1:n ψμ[i] = ms.ψs[i,μ] end

        if !isnothing(F)
            mul!(ms.res,Aμ,ψμ)
            rdiv!(ms.res,ψμ[index])
            @inbounds @simd for i ∈ 1:n F[(μ-1)n+i], F[nm+(μ-1)n+i] = reim(ms.res[i]) end
        end

        if !isnothing(J)
            fill!(J,0)
            ωμϕ = ms.ωs[μ]/ψμ[index]
            ω²μϕ = ms.ωs[μ]^2/ψμ[index]
            @inbounds for ν ∈ 1:m # derivative wrt to ν, evaluated for field μ
                ρνμω²μ = ω²μϕ*ms.ψs[index,ν]
                @fastmath @inbounds @simd for j ∈ 1:n
                    rr,ir = reim(ρνμω²μ*ms.helmholtz.dFχdψr[j,μ,ν]*ψ[j,μ])
                    J[j+(μ-1)n   , j+(ν-1)n] = rr
                    J[j+(μ-1)n+nm, j+(ν-1)n] = ir

                    rr,ir = reim(ρνμω²μ*ms.helmholtz.dFχdψi[j,μ,ν]*ψ[j,μ])
                    J[j+(μ-1)n   , j+(ν-1)n+nm] = rr
                    J[j+(μ-1)n+nm, j+(ν-1)n+nm] = ir
                end

                @fastmath @inbounds @simd for i ∈ 1:n
                    J[(μ-1)n+i,(ν-1)n+index], J[nm+(μ-1)n+i,(ν-1)n+index] = reim(ω²μϕ*(ms.helmholtz.dFχdϕ[i,μ,ν]*ψμ[i]))
                    J[(μ-1)n+i,nm+(ν-1)n+index], J[nm+(μ-1)n+i,nm+(ν-1)n+index] = reim(ω²μϕ*ms.helmholtz.dFχdω[i,μ,ν]*ψμ[i])
                end
                if μ == ν
                    δω = 1e-8
                    ∂∇²ψ∂ω = reduce(+,map((f,σ) -> σ*(f(ms.ωs[μ]+δω)-f(ms.ωs[μ]-δω))*ψμ./(2δω*ψμ[index]),ms.simulation.self_energy.f,ms.helmholtz.Σs[1:end-1]))
                    @fastmath @inbounds for i ∈ 1:n
                        tr, ti = reim(2ωμϕ*(ms.helmholtz.αεpFχ.nzval[i]*ψμ[i]))
                        J[     (μ-1)n+i,nm+(ν-1)n+index] += tr + real(∂∇²ψ∂ω[i])
                        J[nm + (μ-1)n+i,nm+(ν-1)n+index] += ti + imag(∂∇²ψ∂ω[i])
                    end
                end
            end

            rows = rowvals(Aμ)
            vals = nonzeros(Aμ)
            @inbounds for col ∈ 1:n
                if col!==index
                    @inbounds for j ∈ nzrange(Aμ, col)
                        row = rows[j]
                        Ar, Ai = reim(Aμ.nzval[j])
                        J[(μ-1)n+row   , (μ-1)n+col  ] += Ar
                        J[(μ-1)n+nm+row, (μ-1)n+col  ] += Ai
                        J[(μ-1)n+row   , (μ-1)n+nm+col] -= Ai
                        J[(μ-1)n+nm+row, (μ-1)n+nm+col] += Ar
                    end
                end
            end
        end
    end
    return nothing
end

# Helmholtz initialize Jacobian
function initialize_J(ms::HelmholtzSALT)
    n,m = ms.n,ms.m
    nm = n*m
    index = ms.index

    max_count = (16m^2+6m+27m)n
    count = 1
    rows = Vector{Float64}(undef,max_count)
    cols = Vector{Float64}(undef,max_count)
    vals = zeros(Float64,max_count)
    for μ ∈ 1:m
        Aμ = ms.helmholtz(ms.ωs[μ],ms.ωs,ms.ψs)
        for ν ∈ 1:m # derivative wrt to ν, evaluated for field μ
            for j ∈ 1:n
                k = j

                rows[count] = k+(μ-1)n
                cols[count] = j+(ν-1)n          ; count +=1
                rows[count] = k+(μ-1)n+nm
                cols[count] = j+(ν-1)n          ; count +=1

                rows[count] = k+(μ-1)n
                cols[count] = j+(ν-1)n+nm       ; count +=1
                rows[count] = k+(μ-1)n+nm
                cols[count] = j+(ν-1)n+nm       ; count +=1
            end

            for i ∈ 1:n
                rows[count] = (μ-1)n+i
                cols[count] = (ν-1)n+index      ; count +=1

                rows[count] = nm+(μ-1)n+i
                cols[count] = (ν-1)n+index      ; count +=1

                rows[count] = (μ-1)n+i
                cols[count] = nm+(ν-1)n+index   ; count +=1

                rows[count] = nm+(μ-1)n+i
                cols[count] = nm+(ν-1)n+index   ; count +=1
            end
            if μ == ν
                for i ∈ 1:n
                    rows[count] = (μ-1)n+i
                    cols[count] = nm+(ν-1)n+index ; count +=1

                    rows[count] = nm+(μ-1)n+i
                    cols[count] = nm+(ν-1)n+index ; count +=1
                end
            end
        end

        rowsA = rowvals(Aμ)
        for col ∈ 1:n
            if col!==index
                for j ∈ nzrange(Aμ, col)
                    row = rowsA[j]

                    rows[count] = (μ-1)n+row
                    cols[count] = (μ-1)n+col    ; count +=1

                    rows[count] = (μ-1)n+nm+row
                    cols[count] = (μ-1)n+col    ; count +=1

                    rows[count] = (μ-1)n+row
                    cols[count] = (μ-1)n+nm+col ; count +=1

                    rows[count] = (μ-1)n+nm+row
                    cols[count] = (μ-1)n+nm+col ; count +=1
                end
            end
        end
    end
    c = count - 1
    return sparse(rows[1:c],cols[1:c],vals[1:c],2nm,2nm)
end

################################################################################
# Maxwell nonlinear objective + jacobian

# Maxwell
@inline function (ms::MaxwellSALT)(F,J,x::Vector)
    x_to_ωψ!(ms,x)
    n,m = ms.n,ms.m
    nm = n*m
    index = ms.index

    if !isnothing(J) ψx, ψy, ψz = ms.ψs.x, ms.ψs.y, ms.ψs.z end

    ψμ = ms.ψs[:,1]
    @inbounds for μ ∈ 1:m
        Aμ = ms.maxwell(ms.ωs[μ],ms.ωs,ms.ψs)
        @inbounds for i ∈ 1:n ψμ[i] = ms.ψs[i,μ] end

        if !isnothing(F)
            mul!(ms.res,Aμ,ψμ)
            rdiv!(ms.res,ψμ[index])
            @inbounds @simd for i ∈ 1:n F[(μ-1)n+i], F[nm+(μ-1)n+i] = reim(ms.res[i]) end
        end

        if !isnothing(J)
            fill!(J,0)
            ωμϕ = ms.ωs[μ]/ψμ[index]
            ω²μϕ = ms.ωs[μ]^2/ψμ[index]
            @inbounds for ν ∈ 1:m # derivative wrt to ν, evaluated for field μ
                ρνμω²μ = ω²μϕ*ms.ψs[index,ν]
                @fastmath @inbounds @simd for j ∈ 1:n
                    k = mod1(j,n÷3)

                    rr,ir = reim(-ρνμω²μ*ms.maxwell.dFχdψr[j,μ,ν]*ψx[k,μ])
                    J[k+0n÷3+(μ-1)n   , j+(ν-1)n] = rr
                    J[k+0n÷3+(μ-1)n+nm, j+(ν-1)n] = ir

                    rr,ir = reim(-ρνμω²μ*ms.maxwell.dFχdψr[j,μ,ν]*ψy[k,μ])
                    J[k+1n÷3+(μ-1)n   , j+(ν-1)n] = rr
                    J[k+1n÷3+(μ-1)n+nm, j+(ν-1)n] = ir

                    rr,ir = reim(-ρνμω²μ*ms.maxwell.dFχdψr[j,μ,ν]*ψz[k,μ])
                    J[k+2n÷3+(μ-1)n   , j+(ν-1)n] = rr
                    J[k+2n÷3+(μ-1)n+nm, j+(ν-1)n] = ir

                    rr,ir = reim(-ρνμω²μ*ms.maxwell.dFχdψi[j,μ,ν]*ψx[k,μ])
                    J[k+0n÷3+(μ-1)n   , j+(ν-1)n+nm] = rr
                    J[k+0n÷3+(μ-1)n+nm, j+(ν-1)n+nm] = ir

                    rr,ir = reim(-ρνμω²μ*ms.maxwell.dFχdψi[j,μ,ν]*ψy[k,μ])
                    J[k+1n÷3+(μ-1)n   , j+(ν-1)n+nm] = rr
                    J[k+1n÷3+(μ-1)n+nm, j+(ν-1)n+nm] = ir

                    rr,ir = reim(-ρνμω²μ*ms.maxwell.dFχdψi[j,μ,ν]*ψz[k,μ])
                    J[k+2n÷3+(μ-1)n   , j+(ν-1)n+nm] = rr
                    J[k+2n÷3+(μ-1)n+nm, j+(ν-1)n+nm] = ir
                end

                @fastmath @inbounds @simd for i ∈ 1:n
                    J[(μ-1)n+i,(ν-1)n+index], J[nm+(μ-1)n+i,(ν-1)n+index] = reim(-ω²μϕ*(ms.maxwell.dFχdϕ[i,μ,ν]*ψμ[i]))
                    J[(μ-1)n+i,nm+(ν-1)n+index], J[nm+(μ-1)n+i,nm+(ν-1)n+index] = reim(-ω²μϕ*ms.maxwell.dFχdω[i,μ,ν]*ψμ[i])
                end
                if μ == ν
                    @fastmath @inbounds for i ∈ 1:n
                        tr, ti = reim(-2ωμϕ*(ms.maxwell.αεpFχ.nzval[i]*ψμ[i]))
                        J[     (μ-1)n+i,nm+(ν-1)n+index] += tr
                        J[nm + (μ-1)n+i,nm+(ν-1)n+index] += ti
                    end
                end
            end

            rows = rowvals(Aμ)
            vals = nonzeros(Aμ)
            @inbounds for col ∈ 1:n
                if col!==index
                    @inbounds for j ∈ nzrange(Aμ, col)
                        row = rows[j]
                        Ar, Ai = reim(Aμ.nzval[j])
                        J[(μ-1)n+row   , (μ-1)n+col  ] += Ar
                        J[(μ-1)n+nm+row, (μ-1)n+col  ] += Ai
                        J[(μ-1)n+row   , (μ-1)n+nm+col] -= Ai
                        J[(μ-1)n+nm+row, (μ-1)n+nm+col] += Ar
                    end
                end
            end
        end
    end
    return nothing
end

# Maxwell initialize Jacobian
function initialize_J(ms::MaxwellSALT)
    n,m = ms.n,ms.m
    nm = n*m
    index = ms.index

    max_count = (16m^2+6m+27m)n
    count = 1
    rows = Vector{Float64}(undef,max_count)
    cols = Vector{Float64}(undef,max_count)
    vals = zeros(Float64,max_count)
    for μ ∈ 1:m
        Aμ = ms.maxwell(ms.ωs[μ],ms.ωs,ms.ψs)
        for ν ∈ 1:m # derivative wrt to ν, evaluated for field μ
            for j ∈ 1:n
                k = mod1(j,n÷3)

                rows[count] = k+0n÷3+(μ-1)n
                cols[count] = j+(ν-1)n          ; count +=1
                rows[count] = k+0n÷3+(μ-1)n+nm
                cols[count] = j+(ν-1)n          ; count +=1

                rows[count] = k+1n÷3+(μ-1)n
                cols[count] = j+(ν-1)n          ; count +=1
                rows[count] = k+1n÷3+(μ-1)n+nm
                cols[count] = j+(ν-1)n          ; count +=1

                rows[count] = k+2n÷3+(μ-1)n
                cols[count] = j+(ν-1)n          ; count +=1
                rows[count] = k+2n÷3+(μ-1)n+nm
                cols[count] = j+(ν-1)n          ; count +=1

                rows[count] = k+0n÷3+(μ-1)n
                cols[count] = j+(ν-1)n+nm       ; count +=1
                rows[count] = k+0n÷3+(μ-1)n+nm
                cols[count] = j+(ν-1)n+nm       ; count +=1

                rows[count] = k+1n÷3+(μ-1)n
                cols[count] = j+(ν-1)n+nm       ; count +=1
                rows[count] = k+1n÷3+(μ-1)n+nm
                cols[count] = j+(ν-1)n+nm       ; count +=1

                rows[count] = k+2n÷3+(μ-1)n
                cols[count] = j+(ν-1)n+nm       ; count +=1
                rows[count] = k+2n÷3+(μ-1)n+nm
                cols[count] = j+(ν-1)n+nm       ; count +=1
            end

            for i ∈ 1:n
                rows[count] = (μ-1)n+i
                cols[count] = (ν-1)n+index      ; count +=1

                rows[count] = nm+(μ-1)n+i
                cols[count] = (ν-1)n+index      ; count +=1

                rows[count] = (μ-1)n+i
                cols[count] = nm+(ν-1)n+index   ; count +=1

                rows[count] = nm+(μ-1)n+i
                cols[count] = nm+(ν-1)n+index   ; count +=1
            end
            if μ == ν
                for i ∈ 1:n
                    rows[count] = (μ-1)n+i
                    cols[count] = nm+(ν-1)n+index ; count +=1

                    rows[count] = nm+(μ-1)n+i
                    cols[count] = nm+(ν-1)n+index ; count +=1
                end
            end
        end

        rowsA = rowvals(Aμ)
        for col ∈ 1:n
            if col!==index
                for j ∈ nzrange(Aμ, col)
                    row = rowsA[j]

                    rows[count] = (μ-1)n+row
                    cols[count] = (μ-1)n+col    ; count +=1

                    rows[count] = (μ-1)n+nm+row
                    cols[count] = (μ-1)n+col    ; count +=1

                    rows[count] = (μ-1)n+row
                    cols[count] = (μ-1)n+nm+col ; count +=1

                    rows[count] = (μ-1)n+nm+row
                    cols[count] = (μ-1)n+nm+col ; count +=1
                end
            end
        end
    end
    c = count - 1
    return sparse(rows[1:c],cols[1:c],vals[1:c],2nm,2nm)
end

################################################################################
# dictionary between real data vector and fields/frequencies

function ωψ_to_x(ωs::Vector,ψs::VectorField)
    x = Vector{Float64}(undef,2prod(size(ψs)))
    ωψ_to_x!(x,ωs,ψs)
    return x
end

ωψ_to_x!(ms::SALTProblem) = ωψ_to_x!(ms.x,ms.ωs,ms.ψs)
@inline function ωψ_to_x!(x::Vector,ωs::Vector,ψs::VectorField)
    n,m = size(ψs)
    nm = n*m
    index = n÷2+INDEX_OFFSET
    @inbounds for μ ∈ 1:m
        @inbounds for i ∈ 1:n
            if i == index
                x[(μ-1)n+index] = abs(ψs[index,μ])
                x[nm+(μ-1)n+index] = ωs[μ]
            else
                iszero(ψs[index,μ]) ? throw("field vanishes at $index, perhaps due to symmetry") : nothing
                x[(μ-1)n+i], x[nm+(μ-1)n+i] = reim(ψs[i,μ]/ψs[index,μ])
            end
        end
    end
    return nothing
end

x_to_ωψ!(ms::SALTProblem) = x_to_ωψ!(ms.ωs,ms.ψs,ms.x)
x_to_ωψ!(ms::SALTProblem,x) = x_to_ωψ!(ms.ωs,ms.ψs,x)
@inline function x_to_ωψ!(ωs::Vector,ψs::VectorField,x::Vector)
    n,m = size(ψs)
    nm = n*m
    index = n÷2+INDEX_OFFSET
    @inbounds for μ ∈ 1:m
        @inbounds for i ∈ 1:n ψs[i,μ] = complex(x[(μ-1)n+i],x[nm+(μ-1)n+i])*x[(μ-1)n+index] end
        ψs[index,μ] = complex(x[(μ-1)n+index],0)
        ωs[μ] = x[nm+(μ-1)n+index]
    end
    return nothing
end

################################################################################
# the salt bootstrap *UNDER CONSTRUCTION*

struct ThresholdFinder{TSIM,TLS,TKW}
    sim::TSIM
    tls::TLS
    nev::Int
    maxit::Int
    kwargs::TKW
end

function (tf::ThresholdFinder)(D₀::Real)
    tf.tls.D₀ = D₀
    nep = maxwell_nep(tf.sim)
    ω, _ = maxwell_eigen(nep, tf.tls.ωa; nev=tf.nev, maxit=tf.maxit, tf.kwargs...)
    return maximum(imag(ω))
end

function SALTbootstrap(
    sim::Simulation{N},
    tls::TwoLevelSystem;
    nev::Integer=5,
    maxit::Integer=50,
    δ::Real=.1,
    D0::Real=0,
    xatol::Real=1e-6,
    xrtol::Real=0,
    verbose::Bool=false,
    warn::Bool=true,
    iterations::Integer=150,
    kwargs...
    ) where N

    ωs = Vector{Float64}[]
    Ds = Float64[]
    ψs = ElectricField{N}[]

    ωᵗʰ, ψᵗʰ = findthreshold!(tls,sim,D0,δ,warn,verbose,xatol,xrtol,nev,maxit; kwargs...)
    push!(ωs,[ωᵗʰ])
    push!(ψs,0ψᵗʰ)
    push!(Ds,tls.D₀)

    if verbose
        print(stdout,"solving for ")
        printstyled(stdout,"D₀=",fmt("3.3f",tls.D₀),color=PRINTED_COLOR_NUMBER)
        println(stdout)
    end
    tls.D₀ += δ
    salt = maxwell_salt(sim,1)
    nlsolve(salt, ωs[1], δ*ψᵗʰ; iterations=iterations)

    nep = maxwell_nep(sim,salt.ωs,salt.ψs)
    ω, ψ = maxwell_eigen(nep,tls.ωa; nev=nev, maxit=maxit)
    imω,ind = findmax(imag(ω))

    if imω < xatol
        push!(ωs,salt.ωs)
        push!(ψs,salt.ψs)
        push!(Ds,tls.D₀)
        nmodes = 1
    else
        nmodes = 2
    end
    salt = maxwell_salt(sim,nmodes)
    salt.ωs[:] = ω[1]
    salt.ψs[:]
    nlsolve(salt; iterations=iterations)

    nep = maxwell_nep(sim,salt.ωs,salt.ψs)
    ω, ψ = maxwell_eigen(nep,tls.ωa; nev=nev, maxit=maxit)
    imω = maximum(imag(ω))


    return Ds,ωs,ψs
end

function findthreshold!(
            tls::TwoLevelSystem,
            sim::Simulation,
            D0::Real,
            δ::Real,
            warn::Bool,
            verbose::Bool,
            xatol::Real,
            xrtol::Real,
            nev::Integer,
            maxit::Integer;
            kwargs...
            )

    verbose ? println(stdout,"Computing threshold") : nothing
    tls.D₀ = D0
    nep = maxwell_nep(sim)
    ω, _ = maxwell_eigen(nep,tls.ωa; nev=nev, maxit=maxit)
    Δ = abs(-(extrema(real(ω))...))
    if warn && Δ<tls.γp
        printstyled(stdout,"\twarning",color=PRINTED_COLOR_WARN,bold=true)
        printstyled(stdout," frequency spread ", fmt("3.2f",Δ), " < ", fmt("3.2f",2*tls.γp)," = 2γ⟂,",color=PRINTED_COLOR_WARN)
        printstyled(stdout," consider increasing ",color=PRINTED_COLOR_WARN)
        printstyled(stdout,"nev\n",color=PRINTED_COLOR_VARIABLE)
    end
    imω1 = maximum(imag(ω)); imω2 = imω1
    iteration = 1
    verbose ? println(stdout,"\trough estimation") : nothing
    thresh = ThresholdFinder(sim,tls,nev,maxit,kwargs)
    while imω2 < 0
        D0 += δ
        imω1 = imω2
        if verbose
            print(stdout,"\t\titeration ",iteration)
            printstyled(stdout," D₀=",fmt("3.3f",tls.D₀),color=PRINTED_COLOR_NUMBER)
            println(stdout)
            iteration += 1
        end
        imω2 = thresh(D0)
    end
    verbose ? println(stdout,"\trefining threshold") : nothing
    tls.D₀ = find_zero(thresh,(tls.D₀-δ,tls.D₀); xatol=xatol,xrtol=xrtol)
    nep = maxwell_nep(sim)
    ω, ψ = maxwell_eigen(nep,tls.ωa; nev=nev, maxit=maxit, kwargs...)
    _, ind = findmax(imag(ω))
    return real(ω[ind]), ψ(ind)
end

################################################################################
# Pretty Printing

import ..Common.PRINTED_COLOR_LIGHT
import ..Common.PRINTED_COLOR_WARN
import ..Common.PRINTED_COLOR_VARIABLE
import ..Common.PRINTED_COLOR_GOOD
import ..Common.PRINTED_COLOR_BAD
import ..Common.PRINTED_COLOR_INSTRUCTION

function Base.show(io::IO,ms::SALTProblem)
    if !get(io, :SCPA, false)
        print(io,ms.m, " mode")
        ms.m>1 ? print(io,"s") : nothing
        if typeof(ms)<:MaxwellSALT
            printstyled(io," MaxwellSALT",color = PRINTED_COLOR_LIGHT)
        else
            printstyled(io," HelmholtzSALT",color = PRINTED_COLOR_LIGHT)
        end
    end
    if !ms.solved[]
        printstyled(io," (",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"unsolved",color=PRINTED_COLOR_WARN)
        printstyled(io,", pass to ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"NLsolve.nlsolve",color=PRINTED_COLOR_VARIABLE)
        printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
    elseif ms.converged[]
        printstyled(io," (",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"solved",color=PRINTED_COLOR_GOOD)
        printstyled(io,", solution in fields ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"ωs",color=:cyan)
        printstyled(io,", ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"ψs",color=:cyan)
        printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
    else
        printstyled(io," (",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"solve attempted",color=PRINTED_COLOR_BAD)
        printstyled(io,", ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"no convergence",color=PRINTED_COLOR_BAD)
        printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
    end
end

################################################################################
# Plotting

@recipe f(ms::SALTProblem; by=abs2) = ms, by
@recipe f(by::Function,ms::SALTProblem) = ms, by
@recipe f(ms::SALTProblem, by::Function) = ms.ψs, by

@recipe f(sim::Simulation, ms::SALTProblem; by=abs2) = sim, ms, by
@recipe f(ms::SALTProblem, sim::Simulation; by=abs2) = sim, ms, by

@recipe f(ms::SALTProblem, sim::Simulation, by::Function) = sim, ms, by
@recipe f(sim::Simulation, by::Function, ms::SALTProblem) = sim, ms, by
@recipe f(ms::SALTProblem, by::Function, sim::Simulation) = sim, ms, by
@recipe f(by::Function, sim::Simulation, ms::SALTProblem) = sim, ms, by
@recipe f(by::Function, ms::SALTProblem, sim::Simulation) = sim, ms, by

@recipe f(sim::Simulation, ms::SALTProblem, by::Function) = sim, ms.ψs, by

end # module

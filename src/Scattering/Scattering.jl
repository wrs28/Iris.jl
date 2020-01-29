"""
    module Scattering
"""
module Scattering

export HelmholtzLS
export HelmholtzNLS
export MaxwellLS
export MaxwellNLS
export scattering!

files = (
    "1D/Scattering1D.jl",
    # "2D/Scattering2D.jl",
    # "3D/Scattering3D.jl",
    )

interfaces = (
    # "Interfaces/IterativeSolvers.jl",
    )

using ..Common
using ..Spectral
using LinearAlgebra
using NLsolve
using SparseArrays
using RecipesBase

import ..Common.AbstractComplexBL
import ..Common.EQUIVALENT_SOURCE_RELATIVE_CUTOFF
import ..Common.DEFAULT_LUPACK

# include("SPA_Scattering.jl")
# using .SPA_Scattering

################################################################################
# DEFINE STRUCTURES

mutable struct ScatteringSolution{N,M} # general container for scattering solutions
    total::VectorField{N,M}
    incident::VectorField{N,M}
    scattered::VectorField{N,M}
    ω::Float64
    a::Vector{ComplexF64}

    ScatteringSolution(t::T,i::T,s::T,ω::Number,a::AbstractVector) where T<:VectorField{N,M} where {N,M} = new{N,M}(t,i,s,ω,a)
end
# convenience constructor
ScatteringSolution{M}(sim::Simulation,ω::Number,a) where M = ScatteringSolution(VectorField{M}(sim,1),VectorField{M}(sim,1),VectorField{M}(sim,1),ω,a)


mutable struct EquivalentSource{N,M,TW,TA,TSI,TSO,TSIM} # container for source, see dimensional files for details
    incoming_mask::TSI
    outgoing_mask::TSO
    field::VectorField{N,M}
    simulation::TSIM
    ω::TW
    a::TA
    channelflux::Vector{Float64}

    function EquivalentSource(
                incoming_mask::AbstractShape{N},
                outgoing_mask::AbstractShape{N},
                field::VectorField{N,M},
                sim::TSIM,
                ω::TW,
                a::TA,
                channelflux::Vector) where {N,M,TSIM<:Simulation{N},TA,TW}
        return new{N,M,TW,TA,typeof(incoming_mask),typeof(outgoing_mask),TSIM}(
                                    incoming_mask, outgoing_mask,field,sim,ω,a,channelflux)
    end
end

abstract type AbstractScatteringProblem{N,M} end

struct LinearScatteringProblem{N,M,TM,TS} <: AbstractScatteringProblem{N,M}
    operator::TM
    equivalent_source::TS
    j::SparseMatrixCSC{ComplexF64,Int}
    solved::Ref{Bool}
    converged::Ref{Bool}
    solution::ScatteringSolution{N,M}
end

struct NonLinearScatteringProblem{N,M,TLS,TLU} <: AbstractScatteringProblem{N,M} # nonlinear scattering container
    linearscatter::TLS
    ψ::VectorField{N,M}
    residual::Matrix{ComplexF64}
    x::Vector{Float64}
    fixedpoint::Ref{Bool}
    lupack::TLU
end

HelmholtzLS{N} = LinearScatteringProblem{N,1}

HelmholtzNLS{N} = NonLinearScatteringProblem{N,1}

MaxwellLS{N} = LinearScatteringProblem{N,3}

MaxwellNLS{N} = NonLinearScatteringProblem{N,3}


################################################################################
# linear scattering

LinearAlgebra.lu(ls::LinearScatteringProblem, lupack::AbstractLUPACK=DEFAULT_LUPACK) = lu(ls.operator.A,lupack)

"""
    scattering!(ls, [lupack=$(DEFAULT_LUPACK)])

Solve linear scattering problem `ls` with LU package `lupack`.
"""
scattering!(mls::LinearScatteringProblem, lupack::AbstractLUPACK=DEFAULT_LUPACK) = scattering!(mls,lu(mls,lupack))

function scattering!(mls::LinearScatteringProblem, alu::Common.LUfact{TS}) where TS
    if TS<:PSolver
        ldiv!(mls.solution.total.val, alu, Matrix(reshape(j,size(j,1),1)))
    else
        ldiv!(mls.solution.total.val, alu, mls.j)
    end
    for i ∈ eachindex(mls.equivalent_source.field)
        # k = mod1(i,length(mls.j)÷3)
        mls.solution.incident[i] = mls.equivalent_source.field[i]#*mls.equivalent_source.incoming_mask(mls.equivalent_source.field.pos[k])
        mls.solution.scattered[i] = mls.solution.total[i] - mls.solution.incident[i]
    end
    mls.solved[] = true
    mls.converged[] = true
    return nothing
end

foreach(include,interfaces)

################################################################################
# nonlinear scattering NLsolve extension
"""
    scattering!(nls, [init; kwargs...)

Solve nonlinear scattering problem `nls` with (optional) initial guess `init`.
See `NLsolve.nlsolve` for keyword arguments.
"""
scattering!(nls::NonLinearScatteringProblem, args...; kwargs...) = fixedpoint(nls, args...; m=1, kwargs...)

function NLsolve.fixedpoint(nls::NonLinearScatteringProblem{N}, init::VectorField{N}=nls.ψ; kwargs...) where N
    nls.fixedpoint[] = true
    return nlsolve(nls, init; m=1, kwargs..., method=:anderson)
end

function NLsolve.nlsolve(nls::NonLinearScatteringProblem{N,M}; kwargs...) where {N,M}
    if M==1
        mls = HelmholtzLS(nls.sim,nls.ω,nls.a)
    elseif M==3
        mls = MaxwellLS(nls.sim,nls.ω,nls.a)
    end
    scattering!(mls, nls.lupack)
    return nlsolve(nls, mls.solution.total; m=1, kwargs...)
end

function NLsolve.nlsolve(nls::NonLinearScatteringProblem{N}, ψ_init::VectorField{N}; kwargs...) where N
    nls.fixedpoint[] = (get(kwargs,:method,:trust_region)==:anderson) ? true : false
    x_init = ψ_to_x(ψ_init)
    F0 = similar(x_init)
    J0 = nls.fixedpoint[] ? spzeros(Float64,length(F0),length(F0)) : initialize_J(nls)
    df = OnceDifferentiable(nls,nls,nls,x_init,F0,J0)
    results = nlsolve(df, x_init; kwargs...)
    x_to_ψ!(nls,results.zero)
    # @inbounds
    for i ∈ eachindex(nls.ψ)
        # k = mod1(i,length(F0)÷6)
        nls.solution.total[i] = nls.ψ[i]
        nls.solution.incident[i] = nls.equivalent_source.field[i]#*nls.equivalent_source.mask(nls.equivalent_source.field.pos[k])
        nls.solution.scattered[i] = nls.solution.total[i] - nls.solution.incident[i]
    end
    nls.solved[] = true
    nls.converged[] = results.f_converged || results.x_converged
    return results
end

################################################################################
# Nonlinear Scattering Objective and Jacobian

(nls::NonLinearScatteringProblem)(F::Vector,x::Vector) = nls(F,nothing,x)
(nls::NonLinearScatteringProblem)(J::AbstractMatrix,x::Vector) = nls(nothing,J,x)


################################################################################
# Nonlinear Scattering Objective and Jacobian
# Helmholtz Objective and Jacobian

@inline function (nls::HelmholtzNLS)(F,J,x::Vector)
    x_to_ψ!(nls,x)
    if nls.fixedpoint[]
        if !isnothing(F)
            nls.helmholtz(nls.ω,[nls.ω],nls.ψ)
            scattering!(nls.linearscatter,nls.lupack)
            # @fastmath @inbounds @simd
            for i ∈ eachindex(nls.linearscatter.solution.total) F[i], F[length(x)÷2 + i] = reim(nls.solution.total[i] - nls.ψ[i]) end
        end
    else
        A = nls.helmholtz(nls.ω,[nls.ω],nls.ψ)
        if !isnothing(F)
            mul!(nls.residual, A, nls.ψ.val)
            # @fastmath @inbounds @simd
            for i ∈ eachindex(nls.residual) F[i], F[length(x)÷2 + i] = reim(nls.residual[i] - nls.j[i]) end
        end

        if !isnothing(J)
            fill!(J,0)

            N = length(x)÷2
            n = size(nls.ψ,1)
            ψ = nls.ψ
            ω² = nls.ω^2
            # @fastmath
            # @inbounds @simd
            for j ∈ 1:size(nls.ψ,1)
                # k = mod1(j,size(nls.ψ,1)÷3)

                rr,ir = ω².*reim(nls.maxwell.dFχdψr[j,1,1]*ψ[j,1])
                J[j    , j    ] = rr
                J[j + N, j    ] = ir

                rr,ir = ω².*reim(nls.maxwell.dFχdψi[j,1,1]*ψ[j,1])
                J[j    , j + N] = rr
                J[j + N, j + N] = ir
            end

            rows = rowvals(A)
            vals = nonzeros(A)
            for col ∈ 1:n
                # @inbounds @simd
                for j ∈ nzrange(A, col)
                    row = rows[j]
                    Ar, Ai = reim(A.nzval[j])
                    J[row  , col  ] += Ar
                    J[row+N, col  ] += Ai
                    J[row  , col+N] -= Ai
                    J[row+N, col+N] += Ar
                end
            end
        end
    end
    return nothing
end

function initialize_J(nls::HelmholtzNLS)
    A = nls.helmholtz(nls.ω,[nls.ω],nls.ψ)
    n = size(nls.ψ,1)
    max_count = 12n + 27n
    count = 1
    rows = Vector{Float64}(undef,max_count)
    cols = Vector{Float64}(undef,max_count)
    vals = Vector{Float64}(undef,max_count)
    for j ∈ 1:size(nls.ψ,1)
        # k = mod1(j,size(nls.ψ,1)÷3)

        rows[count] = j
        cols[count] = j     ;count+=1
        rows[count] = j+n
        cols[count] = j     ;count+=1

        rows[count] = j
        cols[count] = j + n   ;count+=1
        rows[count] = j+n
        cols[count] = j + n   ;count+=1
    end

    rowsA = rowvals(A)
    for col ∈ 1:n
        for j ∈ nzrange(A, col)
            row = rowsA[j]

            rows[count] = row
            cols[count] = col     ; count+=1

            rows[count] = row + n
            cols[count] = col     ; count+=1

            rows[count] = row
            cols[count] = col + n ; count+=1

            rows[count] = row + n
            cols[count] = col + n ; count+=1
        end
    end
    c = count-1
    return sparse(rows[1:c],cols[1:c],vals[1:c],2n,2n)
end


################################################################################
# Nonlinear Scattering Objective and Jacobian
# Maxwell Objective and Jacobian

@inline function (nls::MaxwellNLS)(F,J,x::Vector)
    x_to_ψ!(nls,x)
    if nls.fixedpoint[]
        if !isnothing(F)
            nls.maxwell(nls.ω,[nls.ω],nls.ψ)
            scattering!(nls.linearscatter,nls.lupack)
            # @fastmath @inbounds @simd
            for i ∈ eachindex(nls.linearscatter.solution.total) F[i], F[length(x)÷2 + i] = reim(nls.solution.total[i] - nls.ψ[i]) end
        end
    else
        A = nls.maxwell(nls.ω,[nls.ω],nls.ψ)
        if !isnothing(F)
            mul!(nls.residual, A, nls.ψ.val)
            # @fastmath @inbounds @simd
            for i ∈ eachindex(nls.residual) F[i], F[length(x)÷2 + i] = reim(nls.residual[i] - nls.j[i]) end
        end

        if !isnothing(J)
            fill!(J,0)

            N = length(x)÷2
            n = size(nls.ψ,1)
            ψx,ψy,ψz = nls.ψ.x,nls.ψ.y,nls.ψ.z
            ω² = nls.ω^2
            # @fastmath
            # @inbounds @simd
            for j ∈ 1:size(nls.ψ,1)
                k = mod1(j,size(nls.ψ,1)÷3)

                rr,ir = ω².*reim(nls.maxwell.dFχdψr[j,1,1]*ψx[k,1])
                J[k+0n÷3    , j    ] = -rr
                J[k+0n÷3 + N, j    ] = -ir

                rr,ir = ω².*reim(nls.maxwell.dFχdψr[j,1,1]*ψy[k,1])
                J[k+1n÷3    , j    ] = -rr
                J[k+1n÷3 + N, j    ] = -ir

                rr,ir = ω².*reim(nls.maxwell.dFχdψr[j,1,1]*ψz[k,1])
                J[k+2n÷3    , j    ] = -rr
                J[k+2n÷3 + N, j    ] = -ir

                rr,ir = ω².*reim(nls.maxwell.dFχdψi[j,1,1]*ψx[k,1])
                J[k+0n÷3    , j + N] = -rr
                J[k+0n÷3 + N, j + N] = -ir

                rr,ir = ω².*reim(nls.maxwell.dFχdψi[j,1,1]*ψy[k,1])
                J[k+1n÷3    , j + N] = -rr
                J[k+1n÷3 + N, j + N] = -ir

                rr,ir = ω².*reim(nls.maxwell.dFχdψi[j,1,1]*ψz[k,1])
                J[k+2n÷3    , j + N] = -rr
                J[k+2n÷3 + N, j + N] = -ir
            end

            rows = rowvals(A)
            vals = nonzeros(A)
            for col ∈ 1:n
                # @inbounds @simd
                for j ∈ nzrange(A, col)
                    row = rows[j]
                    Ar, Ai = reim(A.nzval[j])
                    J[row  , col  ] += Ar
                    J[row+N, col  ] += Ai
                    J[row  , col+N] -= Ai
                    J[row+N, col+N] += Ar
                end
            end
        end
    end
    return nothing
end

function initialize_J(nls::MaxwellNLS)
    A = nls.maxwell(nls.ω,[nls.ω],nls.ψ)
    n = size(nls.ψ,1)
    max_count = 12n + 27n
    count = 1
    rows = Vector{Float64}(undef,max_count)
    cols = Vector{Float64}(undef,max_count)
    vals = Vector{Float64}(undef,max_count)
    for j ∈ 1:size(nls.ψ,1)
        k = mod1(j,size(nls.ψ,1)÷3)

        rows[count] = k+0n÷3
        cols[count] = j     ;count+=1
        rows[count] = k+0n÷3+n
        cols[count] = j     ;count+=1

        rows[count] = k+1n÷3
        cols[count] = j     ;count+=1
        rows[count] = k+1n÷3+n
        cols[count] = j     ;count+=1

        rows[count] = k+2n÷3
        cols[count] = j     ;count+=1
        rows[count] = k+2n÷3+n
        cols[count] = j     ;count+=1

        rows[count] = k+0n÷3
        cols[count] = j + n   ;count+=1
        rows[count] = k+0n÷3+n
        cols[count] = j + n   ;count+=1

        rows[count] = k+1n÷3
        cols[count] = j + n   ;count+=1
        rows[count] = k+1n÷3+n
        cols[count] = j + n   ;count+=1

        rows[count] = k+2n÷3
        cols[count] = j + n   ;count+=1
        rows[count] = k+2n÷3+n
        cols[count] = j + n   ;count+=1
    end

    rowsA = rowvals(A)
    for col ∈ 1:n
        for j ∈ nzrange(A, col)
            row = rowsA[j]

            rows[count] = row
            cols[count] = col     ; count+=1

            rows[count] = row + n
            cols[count] = col     ; count+=1

            rows[count] = row
            cols[count] = col + n ; count+=1

            rows[count] = row + n
            cols[count] = col + n ; count+=1
        end
    end
    c = count-1
    return sparse(rows[1:c],cols[1:c],vals[1:c],2n,2n)
end


################################################################################
# PROPERTY DEFINITIONS

function Base.getproperty(sct::ScatteringSolution, sym::Symbol)
    if Base.sym_in(sym,(:incoming,:inc,:Incident,:Inc,:In,:IN,:in,:i,:I))
        return getfield(sct,:incident)
    elseif Base.sym_in(sym,(:Total,:tot,:Tot,:t,:T))
        return getfield(sct,:total)
    elseif Base.sym_in(sym,(:Scattered,:scat,:Scat,:sct,:SCT,:Sct,:S,:s,:outgoing,:out,:Out,:OUT,:o,:O))
        return getfield(sct,:scattered)
    elseif Base.sym_in(sym,(:frequency,:freq,:omega,:Omega,:Freq,:Frequency,:f,:F))
        return getfield(sct,:ω)
    elseif Base.sym_in(sym,(:amplitude,:amp,:a))
        return getfield(sct,:a)
    elseif Base.sym_in(sym,(:I,:Intensity,:intensity))
        return abs2(getfield(sct,:a))
    else
        return getfield(sct,sym)
    end
end

function Base.propertynames(::ScatteringSolution, private=false)
    if private
        return fieldnames(ScatteringSolution)
    else
        return (:total,:scattered,:incoming,:ω,:a)
    end
end


function Base.getproperty(es::EquivalentSource, sym::Symbol)
    if sym==:sim
        return getfield(es,:simulation)
    else
        return getfield(es,sym)
    end
end


function Base.getproperty(ls::LinearScatteringProblem, sym::Symbol)
    if Base.sym_in(sym,(:sim,:simulation))
        return getfield(getfield(ls,:equivalent_source),:simulation)
    elseif sym==:ω
        return getfield(getfield(ls,:equivalent_source),:ω)
    elseif sym==:a
        return getfield(getfield(ls,:equivalent_source),:a)
    elseif Base.sym_in(sym,(:M,:H,:maxwell,:helmholtz))
        return getfield(ls,:operator)
    elseif Base.sym_in(sym,(:source,))
        return getfield(ls,:equivalent_source)
    elseif Base.sym_in(sym,(:sol,:ψ))
        return getfield(ls,:solution)
    else
        return getfield(ls,sym)
    end
end

function Base.setproperty!(ls::LinearScatteringProblem, sym::Symbol, x)
    if sym==:a
        return setfield!(ls.equivalent_source, :a, complex(x))
    elseif sym==:ω
        return setfield!(ls.equivalent_source, :ω, float(x))
    else
        return setfield!(ls,sym,x)
    end
end

function Base.propertynames(lsp::LinearScatteringProblem, private=false)
    if private
        return fieldnames(LinearScatteringProblem)
    elseif typeof(lsp)<:HelmholtzLS
        return (:simulation, :ω, :a, :equivalent_source, :helmholtz, :solution, :j)
    elseif typeof(lsp)<:MaxwellLS
        return (:simulation, :ω, :a, :equivalent_source, :maxwell, :solution, :j)
    end
end


function Base.getproperty(nls::NonLinearScatteringProblem,sym::Symbol)
    if sym==:res
        return getfield(nls,:residual)
    elseif Base.sym_in(sym,(:sim,:simulation))
        return getproperty(getfield(nls,:linearscatter),:simulation)
    elseif sym==:ω
        return getproperty(getfield(nls,:linearscatter),:ω)
    elseif sym==:a
        return getproperty(getfield(nls,:linearscatter),:a)
    elseif Base.sym_in(sym,(:M,:maxwell,:H,:helmholtz))
        return getproperty(getfield(nls,:linearscatter),:operator)
    elseif Base.sym_in(sym,(:source,:equivalent_source))
        return getproperty(getfield(nls,:linearscatter),:equivalent_source)
    elseif sym==:j
        return getproperty(getfield(nls,:linearscatter),:j)
    elseif Base.sym_in(sym,(:sol,:solution))
        return getproperty(getfield(nls,:linearscatter),:solution)
    elseif sym==:solved
        return getproperty(getfield(nls,:linearscatter),:solved)
    elseif sym==:converged
        return getproperty(getfield(nls,:linearscatter),:converged)
    else
        return getfield(nls,sym)
    end
end

function Base.propertynames(nls::NonLinearScatteringProblem, private=false)
    if private
        return fieldnames(NonLinearScatteringProblem)
    else
        return propertynames(nls.linearscatter,false)
    end
end


################################################################################
# nonlinear transformation between field and re,im vector for nonlinear scattering

function ψ_to_x(ψ::VectorField)
    x = Vector{Float64}(undef,2prod(size(ψ)))
    ψ_to_x!(x,ψ)
    return x
end
ψ_to_x!(nls::NonLinearScatteringProblem, ψ::VectorField) = ψ_to_x!(nls.x,ψ)
@inline function ψ_to_x!(x,ψ::VectorField)
    n,m = size(ψ)
    nm = n*m
    for μ ∈ 1:m
        # @inbounds @simd
        for i ∈ 1:n
            x[(μ-1)n+i], x[nm+(μ-1)n+i] = reim(ψ[i,μ])
        end
    end
    return nothing
end

# function x_to_ψ(sim::Simulation,x,m)
#     ψ = ElectricField(sim.x,m)
#     x_to_ψ!(ψ,x)
#     return ψ
# end
x_to_ψ!(nls::NonLinearScatteringProblem,x) = x_to_ψ!(nls.ψ,x)
@inline function x_to_ψ!(ψ::VectorField,x)
    n,m = size(ψ)
    nm = n*m
    for μ ∈ 1:m
        # @inbounds @simd
        for i ∈ 1:n
            ψ[i,μ] = complex(x[(μ-1)n+i],x[nm+(μ-1)n+i])
        end
    end
    return nothing
end

################################################################################

foreach(include,files)

################################################################################
# Pretty Printing

import ..Common.PRINTED_COLOR_LIGHT
import ..Common.PRINTED_COLOR_GOOD
import ..Common.PRINTED_COLOR_NUMBER
import ..Common.PRINTED_COLOR_BAD
import ..Common.PRINTED_COLOR_WARN
import ..Common.PRINTED_COLOR_VARIABLE
import ..Common.PRINTED_COLOR_INSTRUCTION


function Base.show(io::IO,sc::EquivalentSource{N}) where N
    print(io,"$(N)D ")
    printstyled(io,"EquivalentSource",color=PRINTED_COLOR_LIGHT)
    printstyled(io," (call with argument ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ω",color=PRINTED_COLOR_NUMBER)
    printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
end


function Base.show(io::IO,sc::ScatteringSolution{N}) where N
    print(io,"$(N)D ")
    printstyled(io,"ScatteringSolution",color=PRINTED_COLOR_LIGHT)
    print(io," @ frequency ")
    printstyled(io,"ω=",sc.ω,color=PRINTED_COLOR_NUMBER)
    print(io,", amplitude ")
    printstyled(IOContext(io,:typeinfo=>Array{ComplexF64}),"a=",sc.a,color=PRINTED_COLOR_NUMBER)
end


function Base.show(io::IO,nls::LinearScatteringProblem)
    if typeof(nls)<:HelmholtzLS
        printstyled(io,"HelmholtzLS ",color = PRINTED_COLOR_LIGHT)
    elseif typeof(nls)<:MaxwellLS
        printstyled(io,"MaxwellLS ",color = PRINTED_COLOR_LIGHT)
    end
    print(io,"linear scattering problem")
    print(io," @ frequency ")
    printstyled(io,"ω=",nls.ω,color=PRINTED_COLOR_NUMBER)
    print(io,", amplitudes ")
    printstyled(IOContext(io,:typeinfo=>Array{ComplexF64}),"a=",nls.a,color=PRINTED_COLOR_NUMBER)
    if !nls.solved[]
        printstyled(io," (",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"unsolved",color=PRINTED_COLOR_WARN)
        printstyled(io,", for use in ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"scattering!",color=PRINTED_COLOR_VARIABLE)
        printstyled(io," or with ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"IterativeSolvers",color=PRINTED_COLOR_VARIABLE)
        printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
    elseif nls.converged[]
        printstyled(io," (",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"solved",color=PRINTED_COLOR_GOOD)
        printstyled(io,", solution in field ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"ψ",color=:cyan)
        printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
    else
        print(io," (")
        printstyled(io,"solve attempted",color=PRINTED_COLOR_BAD)
        print(io,", ")
        printstyled(io,"no convergence",color=PRINTED_COLOR_BAD)
        print(io,")")
    end
end


function Base.show(io::IO,nls::NonLinearScatteringProblem)
    if typeof(nls)<:HelmholtzNLS
        printstyled(io,"HelmholtzNLS ",color = PRINTED_COLOR_LIGHT)
    elseif typeof(nls)<:MaxwellNLS
        printstyled(io,"MaxwellNLS ",color = PRINTED_COLOR_LIGHT)
    end
    print(io,"nonlinear scattering problem")
    print(io," @ frequency ")
    printstyled(io,"ω=",nls.ω,color=PRINTED_COLOR_NUMBER)
    print(io,", amplitude ")
    printstyled(IOContext(io,:typeinfo=>Array{ComplexF64}),"a=",nls.a,color=PRINTED_COLOR_NUMBER)
    if !nls.solved[]
        printstyled(io," (",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"unsolved",color=PRINTED_COLOR_WARN)
        printstyled(io,", for use in ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"nlsolve",color=PRINTED_COLOR_VARIABLE)
        printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
    elseif nls.converged[]
        printstyled(io," (",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"solved",color=PRINTED_COLOR_GOOD)
        printstyled(io,", solution in field ",color=PRINTED_COLOR_INSTRUCTION)
        printstyled(io,"ψ",color=:cyan)
        printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
    else
        print(io," (")
        printstyled(io,"solve attempted",color=PRINTED_COLOR_BAD)
        print(io,", ")
        printstyled(io,"no convergence",color=PRINTED_COLOR_BAD)
        print(io,")")
    end
end



################################################################################
# Plotting

@recipe f(ms::AbstractScatteringProblem; by=abs2) = ms, by
@recipe f(by::Function,ms::AbstractScatteringProblem) = ms, by
@recipe f(ms::AbstractScatteringProblem, by::Function) = ms.simulation, ms.solution, by

@recipe f(sim::Simulation, ms::AbstractScatteringProblem; by=abs2) = ms, by
@recipe f(ms::AbstractScatteringProblem, sim::Simulation; by=abs2) = ms, by

@recipe f(ms::AbstractScatteringProblem, sim::Simulation, by::Function) = ms, by
@recipe f(sim::Simulation, by::Function, ms::AbstractScatteringProblem) = ms, by
@recipe f(ms::AbstractScatteringProblem, by::Function, sim::Simulation) = ms, by
@recipe f(by::Function, sim::Simulation, ms::AbstractScatteringProblem) = ms, by
@recipe f(by::Function, ms::AbstractScatteringProblem, sim::Simulation) = ms, by

@recipe f(sim::Simulation, ms::AbstractScatteringProblem, by::Function) = ms.sim, ms.solution, by


@recipe f(ms::ScatteringSolution; by=abs2) = ms, by
@recipe f(by::Function,ms::ScatteringSolution) = ms, by

@recipe f(sim::Simulation, ms::ScatteringSolution; by=abs2) = sim, ms, by
@recipe f(ms::ScatteringSolution, sim::Simulation; by=abs2) = sim, ms, by

@recipe f(ms::ScatteringSolution, sim::Simulation, by::Function) = sim, ms, by
@recipe f(sim::Simulation, by::Function, ms::ScatteringSolution) = sim, ms, by
@recipe f(ms::ScatteringSolution, by::Function, sim::Simulation) = sim, ms, by
@recipe f(by::Function, sim::Simulation, ms::ScatteringSolution) = sim, ms, by
@recipe f(by::Function, ms::ScatteringSolution, sim::Simulation) = sim, ms, by

include("1D/Plotting1D.jl")

end # module

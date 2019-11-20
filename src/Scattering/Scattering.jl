"""
    module Scattering
"""
module Scattering

export maxwell_ls
export maxwell_nls
export scattering
export scattering!
export scattering_nl

# export scattering_spa

files = (
    "1D/Scattering1D.jl",
    # "2D/Scattering2D.jl",
    # "3D/Scattering3D.jl",
    )

using ..Common
using ..Spectral
using IterativeSolvers
using LinearAlgebra
using NLsolve
using SparseArrays

# include("SPA_Scattering.jl")
# using .SPA_Scattering

const DEFAULT_LUPACK = Common.DEFAULT_LUPACK
const EQUIVALENT_SOURCE_RELATIVE_CUTOFF = Common.EQUIVALENT_SOURCE_RELATIVE_CUTOFF
const PRINTED_COLOR_LIGHT = Common.PRINTED_COLOR_LIGHT
const PRINTED_COLOR_GOOD = Common.PRINTED_COLOR_GOOD
const PRINTED_COLOR_NUMBER = Common.PRINTED_COLOR_NUMBER
const PRINTED_COLOR_BAD = Common.PRINTED_COLOR_BAD
const PRINTED_COLOR_WARN = Common.PRINTED_COLOR_WARN
const PRINTED_COLOR_VARIABLE = Common.PRINTED_COLOR_VARIABLE
const PRINTED_COLOR_INSTRUCTION = Common.PRINTED_COLOR_INSTRUCTION


################################################################################
# DEFINE STRUCTURES

mutable struct ScatteringSolution{N} # general container for scattering solutions
    total::ElectricField{N}
    incident::ElectricField{N}
    scattered::ElectricField{N}
    ω::Float64
    a::Vector{ComplexF64}
    ScatteringSolution(t::T,i::T,s::T,ω::Number,a::AbstractVector) where T<:ElectricField{N} where N = new{N}(t,i,s,ω,a)
end

# convenience constructors
ScatteringSolution(sim::Simulation, args...) = ScatteringSolution(sim.x, args...)
ScatteringSolution(pos::Array{T,1},tot::Array,inc::Array,sct::Array,ω,a) where T<:Point =
    ScatteringSolution(ElectricField(pos,tot),ElectricField(pos,inc),ElectricField(pos,sct),ω,a)
ScatteringSolution(pos::Array{T,1},ω::Number,a) where T<:Point =
    ScatteringSolution(ElectricField(pos),ElectricField(pos),ElectricField(pos),ω,a)


mutable struct EquivalentSource{N,TW,TA,TS,TSIM} # container for source, see dimensional files for details
    mask::TS
    field::ElectricField{N}
    simulation::TSIM
    ω::TW
    a::TA
    EquivalentSource(mask::AbstractShape{N},field,sim::TSIM,ω::TW,a::TA) where {N,TSIM<:Simulation,TA,TW} = new{N,TW,TA,typeof(mask),TSIM}(mask,field,sim,ω,a)
end


struct Maxwell_LS{N,TM,TS} # linear scattering problem container
    maxwell::TM
    equivalent_source::TS
    j::SparseMatrixCSC{ComplexF64,Int}
    solved::Ref{Bool}
    converged::Ref{Bool}
    solution::ScatteringSolution{N}
end
maxwell_ls(args...; kwargs...) = Maxwell_LS(args...; kwargs...)


struct Maxwell_NLS{N,TLS,TLU} # nonlinear scattering container
    linearscatter::TLS
    ψ::ElectricField{N}
    residual::Matrix{ComplexF64}
    x::Vector{Float64}
    fixedpoint::Ref{Bool}
    lupack::TLU
end
maxwell_nls(args...;kwargs...) = Maxwell_NLS(args...; kwargs...)


################################################################################
# linear scattering

LinearAlgebra.lu(mls::Maxwell_LS,lupack::AbstractLUPACK=DEFAULT_LUPACK) = lu(mls.maxwell.A,lupack)


"""
    scattering!(::Maxwell_LS, [lupack=$(DEFAULT_LUPACK)])
"""
scattering!(mls::Maxwell_LS, lupack::AbstractLUPACK=DEFAULT_LUPACK) = scattering!(mls,lu(mls,lupack))
function scattering!(mls::Maxwell_LS, alu::Common.LUfact{TS}) where TS
    if TS<:PSolver
        ldiv!(mls.solution.total.val, alu, Matrix(reshape(j,size(j,1),1)))
    else
        ldiv!(mls.solution.total.val, alu, mls.j)
    end
    for i ∈ eachindex(mls.equivalent_source.field)
        k = mod1(i,length(mls.j)÷3)
        mls.solution.incident[i] = mls.equivalent_source.field[i]*mls.equivalent_source.mask(mls.equivalent_source.field.pos[k])
        mls.solution.scattered[i] = mls.solution.total[i] - mls.solution.incident[i]
    end
    mls.solved[] = true
    mls.converged[] = true
    return nothing
end


fnames = (:cg,:minres,:gmres,:idrs,:bicgstabl,:chebyshev,:jacobi,:gauss_seidel,:sor,:ssor,:lsmr,:lsqr)
for fname ∈ fnames
    @eval begin
        function IterativeSolvers.$(fname)(mls::Maxwell_LS, args...; kwargs...)
            $(Symbol(fname,"!"))(mls.solution.tot,mls.maxwell.A,mls.j, args...; kwargs...)
            return mls.solution
        end
        function IterativeSolvers.$(Symbol(fname,"!"))(mls::Maxwell_LS, args...; kwargs...)
            _, hist = $(Symbol(fname,"!"))(mls.solution.tot,mls.maxwell.A,mls.j, args...; kwargs..., log=true)
            for i ∈ eachindex(mls.equivalent_source.field)
                k = mod1(i,length(j)÷3)
                mls.solution.incident[i] = mls.equivalent_source.field[i]*mls.equivalent_source.mask(mls.equivalent_source.field.pos[k])
                mls.solution.scattered[i] = mls.solution.total[i] - mls.solution.incident[i]
            end
            mls.solved[] = true
            mls.converged[] = hist.isconverged
            return nothing
        end
    end
end


################################################################################
# nonlinear scattering

function NLsolve.fixedpoint(nls::Maxwell_NLS{N}, init::ElectricField{N}=nls.ψ; kwargs...) where N
    nls.fixedpoint[] = true
    return nlsolve(nls, init; method=:anderson, kwargs...)
end

function NLsolve.nlsolve(nls::Maxwell_NLS; kwargs...)
    mls = Maxwell_LS(nls.sim,nls.ω,nls.a)
    scattering!(mls, nls.lupack)
    return nlsolve(nls, mls.solution.tot; kwargs...)
end

function NLsolve.nlsolve(nls::Maxwell_NLS{N}, ψ_init::ElectricField{N}; kwargs...) where N
    nls.fixedpoint[] = (get(kwargs,:method,:trust_region)==:anderson) ? true : false
    x_init = ψ_to_x(ψ_init)
    F0 = similar(x_init)
    J0 = nls.fixedpoint[] ? spzeros(Float64,length(F0),length(F0)) : initialize_J(nls)
    df = OnceDifferentiable(nls,nls,nls,x_init,F0,J0)
    results = nlsolve(df, x_init; kwargs...)
    x_to_ψ!(nls,results.zero)
    @inbounds for i ∈ eachindex(nls.ψ)
        k = mod1(i,length(F0)÷6)
        nls.solution.total[i] = nls.ψ[i]
        nls.solution.incident[i] = nls.equivalent_source.field[i]*nls.equivalent_source.mask(nls.equivalent_source.field.pos[k])
        nls.solution.scattered[i] = nls.solution.total[i] - nls.solution.incident[i]
    end
    nls.solved[] = true
    nls.converged[] = results.f_converged || results.x_converged
    return results
end

(nls::Maxwell_NLS)(F::Vector,x::Vector) = nls(F,nothing,x)
(nls::Maxwell_NLS)(J::AbstractMatrix,x::Vector) = nls(nothing,J,x)
@inline function (nls::Maxwell_NLS)(F,J,x::Vector)
    x_to_ψ!(nls,x)
    if nls.fixedpoint[]
        if !isnothing(F)
            nls.maxwell(nls.ω,[nls.ω],nls.ψ)
            scattering!(nls.linearscatter,nls.lupack)
            @fastmath @inbounds @simd for i ∈ eachindex(nls.linearscatter.solution.total) F[i], F[length(x)÷2 + i] = reim(nls.solution.total[i] - nls.ψ[i]) end
        end
    else
        A = nls.maxwell(nls.ω,[nls.ω],nls.ψ)
        if !isnothing(F)
            mul!(nls.residual, A, nls.ψ.val)
            @fastmath @inbounds @simd for i ∈ eachindex(nls.residual) F[i], F[length(x)÷2 + i] = reim(nls.residual[i] - nls.j[i]) end
        end

        if !isnothing(J)
            fill!(J,0)

            N = length(x)÷2
            n = size(nls.ψ,1)
            ψx,ψy,ψz = nls.ψ.x,nls.ψ.y,nls.ψ.z
            ω² = nls.ω^2
            # @fastmath
            @inbounds @simd for j ∈ 1:size(nls.ψ,1)
                k = mod1(j,size(nls.ψ,1)÷3)

                rr,ir = ω².*reim(nls.M.αdχdψr[j,1,1]*ψx[k,1])
                J[k+0n÷3    , j    ] = -rr
                J[k+0n÷3 + N, j    ] = -ir

                rr,ir = ω².*reim(nls.M.αdχdψr[j,1,1]*ψy[k,1])
                J[k+1n÷3    , j    ] = -rr
                J[k+1n÷3 + N, j    ] = -ir

                rr,ir = ω².*reim(nls.M.αdχdψr[j,1,1]*ψz[k,1])
                J[k+2n÷3    , j    ] = -rr
                J[k+2n÷3 + N, j    ] = -ir

                rr,ir = ω².*reim(nls.M.αdχdψi[j,1,1]*ψx[k,1])
                J[k+0n÷3    , j + N] = -rr
                J[k+0n÷3 + N, j + N] = -ir

                rr,ir = ω².*reim(nls.M.αdχdψi[j,1,1]*ψy[k,1])
                J[k+1n÷3    , j + N] = -rr
                J[k+1n÷3 + N, j + N] = -ir

                rr,ir = ω².*reim(nls.M.αdχdψi[j,1,1]*ψz[k,1])
                J[k+2n÷3    , j + N] = -rr
                J[k+2n÷3 + N, j + N] = -ir
            end

            rows = rowvals(A)
            vals = nonzeros(A)
            for col ∈ 1:n
                @inbounds @simd for j ∈ nzrange(A, col)
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

function initialize_J(nls::Maxwell_NLS)
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

function Base.getproperty(sct::ScatteringSolution,sym::Symbol)
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


function Base.getproperty(es::EquivalentSource,sym::Symbol)
    if sym==:sim
        return getfield(es,:simulation)
    else
        return getfield(es,sym)
    end
end

function Base.propertynames(::EquivalentSource, private=false)
    if private
        return fieldnames(EquivalentSource)
    else
        return (:mask, :field, :simulation, :ω, :a)
    end
end


function Base.getproperty(ls::Maxwell_LS,sym::Symbol)
    if Base.sym_in(sym,(:sim,:simulation))
        return getfield(getfield(ls,:equivalent_source),:simulation)
    elseif sym==:ω
        return getfield(getfield(ls,:equivalent_source),:ω)
    elseif sym==:a
        return getfield(getfield(ls,:equivalent_source),:a)
    elseif Base.sym_in(sym,(:M,))
        return getfield(ls,:maxwell)
    elseif Base.sym_in(sym,(:source,))
        return getfield(ls,:equivalent_source)
    elseif Base.sym_in(sym,(:sol,))
        return getfield(ls,:solution)
    else
        return getfield(ls,sym)
    end
end

function Base.propertynames(::Maxwell_LS, private=false)
    if private
        return fieldnames(Maxwell_LS)
    else
        return (:simulation, :ω, :a, :equivalent_source, :maxwell, :solution, :j)
    end
end


function Base.getproperty(nls::Maxwell_NLS,sym::Symbol)
    if sym==:res
        return getfield(nls,:residual)
    elseif Base.sym_in(sym,(:sim,:simulation))
        return getproperty(getfield(nls,:linearscatter),:simulation)
    elseif sym==:ω
        return getproperty(getfield(nls,:linearscatter),:ω)
    elseif sym==:a
        return getproperty(getfield(nls,:linearscatter),:a)
    elseif Base.sym_in(sym,(:M,:maxwell))
        return getproperty(getfield(nls,:linearscatter),:maxwell)
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

function Base.propertynames(nls::Maxwell_NLS,private=false)
    if private
        return fieldnames(Maxwell_NLS)
    else
        return propertynames(nls.linearscatter,false)
    end
end

################################################################################

foreach(include,files)

################################################################################

function ψ_to_x(ψ::ElectricField)
    x = Vector{Float64}(undef,2prod(size(ψ)))
    ψ_to_x!(x,ψ)
    return x
end
ψ_to_x!(nls::Maxwell_NLS,ψ::ElectricField) = ψ_to_x!(nls.x,ψ)
@inline function ψ_to_x!(x,ψ::ElectricField)
    n,m = size(ψ)
    nm = n*m
    for μ ∈ 1:m
        @inbounds @simd for i ∈ 1:n
            x[(μ-1)n+i], x[nm+(μ-1)n+i] = reim(ψ[i,μ])
        end
    end
    return nothing
end

function x_to_ψ(sim::Simulation,x,m)
    ψ = ElectricField(sim.x,m)
    x_to_ψ!(ψ,x)
    return ψ
end
x_to_ψ!(nls::Maxwell_NLS,x) = x_to_ψ!(nls.ψ,x)
@inline function x_to_ψ!(ψ::ElectricField,x)
    n,m = size(ψ)
    nm = n*m
    for μ ∈ 1:m
        @inbounds @simd for i ∈ 1:n
            ψ[i,μ] = complex(x[(μ-1)n+i],x[nm+(μ-1)n+i])
        end
    end
    return nothing
end

################################################################################
# Pretty Printing

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


function Base.show(io::IO,nls::Maxwell_LS)
    printstyled(io,"Maxwell_LS ",color = PRINTED_COLOR_LIGHT)
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


function Base.show(io::IO,nls::Maxwell_NLS)
    printstyled(io,"Maxwell_NLS ",color = PRINTED_COLOR_LIGHT)
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

end # module

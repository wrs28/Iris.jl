"""
	module Spectral
"""
module Spectral

export maxwell_lep
export maxwell_cf
export maxwell_nep
export maxwell_eigs
export maxwell_eigs_l
export maxwell_eigs_cf
export maxwell_eigs_nl
export IrisLinSolverCreator

files = ("1D/Spectral1D.jl",
        # "2D/Spectral2D.jl",
        # "3D/Spectral3D.jl"
        )

using ..Common
using ArnoldiMethod
using ArnoldiMethodTransformations: AbstractSolver, MSolver, PSolver, USolver, initialize_according_to_package
using Flatten
using LinearAlgebra
using NonlinearEigenproblems
using Pardiso
using RecipesBase
using SparseArrays

const PRINTED_COLOR_DARK = Common.PRINTED_COLOR_DARK
const PRINTED_COLOR_VARIABLE = Common.PRINTED_COLOR_VARIABLE


# linear
################################################################################
struct Maxwell_LEP{N,TM}
    M::TM
    Maxwell_LEP(m::Maxwell{N}) where N = new{N,typeof(m)}(m)
end

Maxwell_LEP(args...;kwargs...) = Maxwell_LEP(Maxwell(args...;kwargs...))

function (lep::Maxwell_LEP)(ω)
    lep.M(ω)
    return lep.M.D², lep.M.αεpχ
end

function maxwell_eigs(
            lep::Maxwell_LEP,
            ω::Number;
            verbose::Bool=false,
            kwargs...)

    A,B = lep(ω)
    decomp, history = partialschur(A, B, ω^2; diag_inv_B=true, kwargs...)
    history.converged || @warn "$(history.nev - history.nconverged) evecs did not converge"
    verbose ? println(history) : nothing

    λ::Vector{ComplexF64}, v::Matrix{ComplexF64} = partialeigen(decomp, ω^2)
    normalize!(lep.M.sim,v,B) # Normalize according to (ψ₁,ψ₂)=δ₁₂
    e = ElectricField(lep.M.sim.x,v)
    return sqrt.(λ), e
end


# CF
################################################################################
struct Maxwell_CF{N,TM}
    M::TM
    Maxwell_CF(m::Maxwell{N}) where N = new{N,typeof(m)}(m)
end

Maxwell_CF(args...;kwargs...) = Maxwell_CF(Maxwell(args...;kwargs...))

function (cf::Maxwell_CF)(ω)
    cf.M(ω)
    ω² = ω^2
    rows = rowvals(cf.M.A)
    vals = nonzeros(cf.M.A)
    _, n = size(cf.M.A)
    for i ∈ 1:n
        col = i
        for j ∈ nzrange(cf.M.A, i)
            row = rows[j]
            vals[j] = cf.M.D²[row,col] - ω²*cf.M.αε[row,col]
        end
    end
    return cf.M.A, spdiagm(0=>repeat(ω²*cf.M.sim.α[1].*cf.M.sim.F,3))
end

function maxwell_eigs(
            cf::Maxwell_CF,
            ω::Number;
            η::Number=0,
            verbose::Bool=false,
            kwargs...)

    A,B = cf(ω)
    decomp, history = partialschur(A, B, η; diag_inv_B=false, kwargs...)
    history.converged || @warn "$(history.nev-history.nconverged) eigenvectors failed to converge"
    verbose ? println(history) : nothing

    λ::Array{ComplexF64,1}, v::Array{ComplexF64,2} = partialeigen(decomp,η)
    normalize!(cf.M.sim,v,B)
    e = ElectricField(cf.M.sim.x,v)
    return λ, e
end


# nonlinear ep
################################################################################
struct Maxwell_NEP{N,TM,TNEP}
    M::TM
    nep::TNEP
    Maxwell_NEP(m::Maxwell{N},NEP) where N = new{N,typeof(m),typeof(NEP)}(m,NEP)
end

(mnep::Maxwell_NEP)(ω) = compute_Mder(mnep.nep,ω)

fnames = names(NonlinearEigenproblems.NEPSolver)
for fn ∈ fnames
    @eval NonlinearEigenproblems.$(fn)(nep::Maxwell_NEP,args...;kwargs...) = $(fn)(nep.nep,args...;kwargs...)
end

function maxwell_eigs(
            nep::Maxwell_NEP,
            ω::Number;
            verbose::Bool = false,
            lupack::AbstractSolver = USolver(),
            logger::Integer = Int(verbose),
            linsolvercreator::LinSolverCreator = IrisLinSolverCreator(lupack),
            kwargs...
            ) where S<:AbstractSolver

    λ::Vector{ComplexF64}, v::Matrix{ComplexF64} = iar(nep; σ=ω,logger=logger,linsolvercreator=linsolvercreator, kwargs...)
    e = ElectricField(sim.x,v)
    return λ, e
end


struct IrisLinSolverCreator{S} <: LinSolverCreator solver::S end
IrisLinSolverCreator(s::T) where T<:AbstractSolver = IrisLinSolverCreator{T}(s)
IrisLinSolverCreator(::Type{T}) where T<:AbstractSolver = IrisLinSolverCreator{T}(T())
IrisLinSolverCreator() where T<:AbstractSolver = IrisLinSolverCreator(USolver)

struct IrisLinSolver{S,L,T} <: LinSolver
    LU::L
    x::Vector{T}
    Z::SparseMatrixCSC{T,Int}
end

function NonlinearEigenproblems.create_linsolver(ILS::IrisLinSolverCreator{S},nep,λ) where S
    M = compute_Mder(nep,λ)
    type = eltype(M)
    x = Vector{type}(undef,size(M,1))
    temp2 = Vector{type}(undef,size(M,1))
    _,LU = initialize_according_to_package(ILS.solver,M,issymmetric(M),type,x,temp2)
    return IrisLinSolver{S,typeof(LU),eltype(x)}(LU,x,spzeros(type,length(x),length(x)))
end

function NonlinearEigenproblems.lin_solve(solver::IrisLinSolver{S},b::Vector;tol=eps()) where S
    if S<:PSolver
        pardiso(solver.LU,solver.x,solver.Z,b)
    else
        ldiv!(solver.x, solver.LU, b)
    end
    return solver.x
end

################################################################################

# include 1D,2D,3D
foreach(include,files)


# Pretty Printing
################################################################################
function Base.show(io::IO,::Maxwell_LEP{N}) where N
    print(io,N,"D ")
    printstyled(io,"Maxwell_LEP ",color=PRINTED_COLOR_DARK)
    print(io,"Linear Eigenproblem (call w/ args ")
    printstyled(io,"ω",color=PRINTED_COLOR_VARIABLE)
    print(io,", [")
    printstyled(io,"ωs",color=PRINTED_COLOR_VARIABLE)
    print(io,",")
    printstyled(io,"ψs",color=PRINTED_COLOR_VARIABLE)
    print(io,"] -> ")
    printstyled(io,"A",color=PRINTED_COLOR_VARIABLE)
    print(io,", ")
    printstyled(io,"B",color=PRINTED_COLOR_VARIABLE)
    print(io,", where ")
    printstyled(io,"Aψ=ω²Bψ",color=PRINTED_COLOR_VARIABLE)
    print(io,")")
end


function Base.show(io::IO,::Maxwell_CF{N}) where N
    print(io,N,"D ")
    printstyled(io,"Maxwell_CF ",color=PRINTED_COLOR_DARK)
    print(io,"Linear Eigenproblem (call w/ args ")
    printstyled(io,"ω",color=PRINTED_COLOR_VARIABLE)
    print(io,", [")
    printstyled(io,"ωs",color=PRINTED_COLOR_VARIABLE)
    print(io,",")
    printstyled(io,"ψs",color=PRINTED_COLOR_VARIABLE)
    print(io,"] -> ")
    printstyled(io,"A",color=PRINTED_COLOR_VARIABLE)
    print(io,", ")
    printstyled(io,"B",color=PRINTED_COLOR_VARIABLE)
    print(io,", where ")
    printstyled(io,"Aψ=ω²Bψ",color=PRINTED_COLOR_VARIABLE)
    print(io,")")
end


function Base.show(io::IO,::Maxwell_NEP{N}) where N
    print(io,N,"D ")
    printstyled(io,"Maxwell_NEP ",color=PRINTED_COLOR_DARK)
    print(io,"Nonlinear Eigenproblem for use with ")
    printstyled(io,"nlsolve",color=PRINTED_COLOR_VARIABLE)
    print(io," (or evaluate operator by calling w/arg ")
    printstyled(io,"ω",color=PRINTED_COLOR_VARIABLE)
    print(io,")")
end


# shortcuts
################################################################################
"""
    maxwell_lep(::Simulation{1}; ky=0, kz=0)
"""
maxwell_lep(args...;kwargs...) = Maxwell_LEP(args...;kwargs...)


"""
    maxwell_lep(::Simulation{1}; ky=0, kz=0)
"""
maxwell_cf(args...;kwargs...) = Maxwell_CF(args...;kwargs...)

"""
    maxwell_nep(::Simulation{1}; ky=0, kz=0, kwargs...)
"""
maxwell_nep(args...;kwargs...) = Maxwell_NEP(args...;kwargs...)


function LinearAlgebra.normalize!(sim::Simulation{1},ψ,B)
    n = length(sim.domain_indices)
    lattices = map(d->d.lattice,sim.domains)
    dx = map(i->lattices[i].dx,sim.domain_indices)
    for i ∈ 1:size(ψ,2)
        𝒩² = zero(eltype(ψ))
        for j ∈ 1:size(ψ,1)
            𝒩² += (ψ[j,i].^2)*B[j,j]*dx[mod1(j,n)]
        end
        ψ[:,i] /= sqrt(𝒩²)*exp(complex(0,angle(ψ[end÷2-1,i])))
    end
    return nothing
end

end # module

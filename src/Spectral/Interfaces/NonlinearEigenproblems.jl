using NonlinearEigenproblems
using LinearAlgebra

export IrisLinSolverCreator

using Pardiso
using ArnoldiMethodTransformations: PSolver, initialize_according_to_package

# overload all the NonlinearEigenproblems solvers
fnames = names(NonlinearEigenproblems.NEPSolver)
for fn ∈ fnames
    @eval begin
        """
            $($fn)(nep::MaxwellNEP, args...; kwargs...) -> λ, e::ElectricField
        """
        function NonlinearEigenproblems.$(fn)(nep::MaxwellNEP, args...; kwargs...)
            λ, v = $(fn)(nep.nep,args...;kwargs...)
            return λ, ElectricField(nep,v)
        end
    end
end

################################################################################
# DEFINE LINSOLVERCREATOR

struct IrisLinSolverCreator{S} <: LinSolverCreator solver::S end
IrisLinSolverCreator(s::T) where T<:AbstractLUPACK = IrisLinSolverCreator{T}(s)
IrisLinSolverCreator(::Type{T}) where T<:AbstractLUPACK = IrisLinSolverCreator{T}(T())
IrisLinSolverCreator() where T<:AbstractLUPACK = IrisLinSolverCreator(DEFAULT_LUPACK)

struct IrisLinSolver{S,L,T} <: LinSolver
    LU::L
    x::Vector{T}
    Z::SparseMatrixCSC{T,Int}
end

function NonlinearEigenproblems.create_linsolver(ILS::IrisLinSolverCreator{S}, nep, λ) where S
    LU = lu(compute_Mder(nep,λ),S)
    x = Vector{type}(undef,size(M,1))
    return IrisLinSolver{S,typeof(LU),eltype(x)}(LU,x,spzeros(type,length(x),length(x)))
end

function NonlinearEigenproblems.lin_solve(solver::IrisLinSolver{S}, b::Vector; tol=eps()) where S
    ldiv!(solver.x, solver.LU, b)
    return solver.x
end

################################################################################
# maxwelleigen

"""
	maxwelleigen_nonlineareigenproblems(nep::MaxwellNEP, ω; [nev=1, verbose=false, lupack=$DEFAULT_LUPACK, linsolvercreator=IrisLinSolverCreator(lupack)]) -> ωs, ψs::ElectricField
"""
function maxwelleigen_nonlineareigenproblems(
            nep::MaxwellNEP,
            ω::Number;
            nev::Int=1,
            verbose::Bool = false,
            lupack::AbstractLUPACK = DEFAULT_LUPACK,
            logger::Integer = Int(verbose),
            linsolvercreator::LinSolverCreator = IrisLinSolverCreator(lupack),
            kwargs...
            ) where S<:AbstractLUPACK

    neigs = get(kwargs,:neigs,nev)
    λ::Vector{ComplexF64}, v::Matrix{ComplexF64} = iar_chebyshev(nep; neigs=neigs, σ=ω,logger=logger, linsolvercreator=linsolvercreator, kwargs...)
    return λ, ElectricField(nep.simulation,v)
end

module NonlinearEigenproblemsInterface

using NonlinearEigenproblems
using LinearAlgebra
using SparseArrays
using ..Common

import ..AbstractEigenproblem
import ..AbstractLinearEigenproblem
import ..AbstractCFEigenproblem
import ..AbstractNonlinearEigenproblem
import ..DEFAULT_NONLINEAR_EIGENSOLVER
import ..INDEX_OFFSET
import ..HelmholtzNEP
import ..MaxwellNEP
import ..DOC_NEP_H
import ..DOC_NEP_M
import ..helmholtzeigen
import ..maxwelleigen
import ..orthogonalize!

export IrisLinSolverCreator

# overload all the NonlinearEigenproblems solvers
fnames = names(NonlinearEigenproblems.NEPSolver)
for fn ∈ fnames
    @eval begin
        """
            $($fn)(nep::HelmholtzNEP, args...; kwargs...) -> λ, s::ScalarField
        """
        function NonlinearEigenproblems.$(fn)(nep::HelmholtzNEP, args...; kwargs...)
            λ, v = $(fn)(nep.nep,args...;kwargs...)
            return λ, ScalarField(nep,v)
        end
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
"""
	IrisLinSolverCreator{S}
"""
struct IrisLinSolverCreator{S} <: LinSolverCreator solver::S end
"""
	IrisLinSolverCreator([lupack])

`lupack` defaults to $DEFAULT_LUPACK.
For use in [`NonlinearEigenproblems`](https://github.com/nep-pack/NonlinearEigenproblems.jl).
"""
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
    x = Vector{ComplexF64}(undef,size(nep,1))
    return IrisLinSolver{S,typeof(LU),eltype(x)}(LU,x,spzeros(eltype(x),length(x),length(x)))
end

function NonlinearEigenproblems.lin_solve(solver::IrisLinSolver{S}, b::Vector; tol=eps()) where S
    ldiv!(solver.x, solver.LU, b)
    return solver.x
end

################################################################################
# maxwelleigen

function iris_eigen_nonlineareigenproblems(
            nep::AbstractNonlinearEigenproblem,
            ω::Number,
            args...;
            verbose::Bool = false,
            lupack::AbstractLUPACK = DEFAULT_LUPACK,
            linsolvercreator::LinSolverCreator = IrisLinSolverCreator(lupack),
            kwargs...
            ) where S<:AbstractLUPACK

    isempty(args) ? nothing : nep(ω,args...)
    λ::Vector{ComplexF64}, v::Matrix{ComplexF64} = iar_chebyshev(nep; σ=ω,logger=Int(verbose), linsolvercreator=linsolvercreator, kwargs...)
	if typeof(nep) <: HelmholtzNEP
		ψ = ScalarField(nep,v)
	elseif typeof(nep) <: MaxwellNEP
		ψ = ElectricField(nep,v)
	end
    normalize!(nep.simulation, ψ, nep.nep.A[end-1], size(ψ,1)÷2+INDEX_OFFSET) # Normalize according to (ψ₁,ψ₂)=δ₁₂
	orthogonalize!(ψ,nep.simulation, λ, nep.nep.A[end-1], nep.ky, nep.kz)
    return λ,ψ
end

if DEFAULT_NONLINEAR_EIGENSOLVER == :NonlinearEigenproblems
	@doc """
	$DOC_NEP_H
	`neigs` Number of eigenpairs (`6`);
	`maxit` Maximum iterations (`30`);
	`verbose` (`false`);
	`v` Initial eigenvector approximation (`rand(size(nep,1),1)`)

	Advanced keywords:
	`logger` Integer print level (`verbose`);
	`tolerance` Convergence tolerance (`eps()*10000`);
	`lupack` LU solver (`$DEFAULT_LUPACK`);
	`linsolvercreator` Define linear solver (`IrisLinSolverCreator(lupack)`);
	`check_error_every` Check for convergence every this number of iterations (`1`);
	`γ` Eigenvalue scaling in `λ=γs+Ω` (`1`);
	`orthmethod` Orthogonalization method, see `IterativeSolvers` (`DGKS`);
	`errmeasure` See `NonlinearEigenproblems` documentation;
	`compute_y0_method` How to find next vector in Krylov space, see `NonlinearEigenproblems` docs (`ComputeY0ChebAuto`)
	a (`-1`);
	b (`1`);
	""" ->
	helmholtzeigen(nep::HelmholtzNEP,args...;kwargs...) = iris_eigen_nonlineareigenproblems(nep,args...;kwargs...)

	@doc """
	$DOC_NEP_M
	`neigs` Number of eigenpairs (`6`);
	`maxit` Maximum iterations (`30`);
	`verbose` (`false`);
	`v` Initial eigenvector approximation (`rand(size(nep,1),1)`)

	Advanced keywords:
	`logger` Integer print level (`verbose`);
	`tolerance` Convergence tolerance (`eps()*10000`);
	`lupack` LU solver (`$DEFAULT_LUPACK`);
	`linsolvercreator` Define linear solver (`IrisLinSolverCreator(lupack)`);
	`check_error_every` Check for convergence every this number of iterations (`1`);
	`γ` Eigenvalue scaling in `λ=γs+Ω` (`1`);
	`orthmethod` Orthogonalization method, see `IterativeSolvers` (`DGKS`);
	`errmeasure` See `NonlinearEigenproblems` documentation;
	`compute_y0_method` How to find next vector in Krylov space, see `NonlinearEigenproblems` docs (`ComputeY0ChebAuto`)
	a (`-1`);
	b (`1`);
	""" ->
	maxwelleigen(nep::MaxwellNEP,args...;kwargs...) = iris_eigen_nonlineareigenproblems(nep,args...;kwargs...)
end

end


using .NonlinearEigenproblemsInterface

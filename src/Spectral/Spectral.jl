"""
	module Spectral

For the computation of various linear and nonlinear eigenvalue problems.

	*LEP: Linear Eigenvalue Problem, i.e. no bulk dispersion, fixed boundary conditions
		(can use with PMLs/cPMLs for resonance or RSM calculations)
	*CF: Constant Flux states, i.e. at a given frequency, find eigen-susceptibilities and
		associated eigenfunctions. Can use either matched boundaries or PMLs. The CF problem
		is a particular kind of linear eigenvalue problem, though here it is treated separately.
	*NEP: Nonlinear Eigenvalue Problem, i.e. with bulk dispersion, open boundary conditions,
		or both. Used to compute resonances/RSMs, no need for PMLs.

	API: `MaxwellLEP`, `MaxwellCF`, `MaxwellNEP`, for use with `maxwelleigen` (high-level),
		or directly with `ArnoldiMethod.jl`, `Arpack.jl`, `IterativeSolvers.jl`,
		`KrylovKit.jl`, `LinearAlgebra`, `NonlinearEigenproblems.jl` (low-level).
"""
module Spectral

export MaxwellLEP
export MaxwellCF
export MaxwellNEP
export maxwelleigen

dimensions = (
		"1D/Spectral1D.jl",
        # "2D/Spectral2D.jl",
        # "3D/Spectral3D.jl"
        )

interfaces = (
	"Interfaces/ArnoldiMethod.jl",
	"Interfaces/Arpack.jl",
	"Interfaces/IterativeSolvers.jl",
	"Interfaces/KrylovKit.jl",
	"Interfaces/LinearAlgebra.jl",
	"Interfaces/NonlinearEigenproblems.jl"
	)

using ..Common
using Flatten

import ..Common.DEFAULT_LINEAR_EIGENSOLVER
import ..Common.DEFAULT_NONLINEAR_EIGENSOLVER
import ..Common.INDEX_OFFSET
import ..Common.PRINTED_COLOR_DARK
import ..Common.PRINTED_COLOR_VARIABLE
import ..Common.PRINTED_COLOR_INSTRUCTION
import ..Common.DEFAULT_LUPACK

abstract type AbstractMaxwellEigenproblem{N} end
abstract type AbstractMaxwellLinearEigenproblem{N} <: AbstractMaxwellEigenproblem{N} end
abstract type AbstractMaxwellNonlinearEigenproblem{N} <: AbstractMaxwellEigenproblem{N} end

# Linear Eigenvalue Problem
struct MaxwellLEP{N,TM} <: AbstractMaxwellLinearEigenproblem{N}
    maxwell::TM
	εpχ::Vector{ComplexF64}
    MaxwellLEP(maxwell::Maxwell{N}) where N = new{N,typeof(maxwell)}(maxwell,maxwell.αεpχ)
end

# Linear Eigenvalue CF Problem
struct MaxwellCF{N,TM} <: AbstractMaxwellLinearEigenproblem{N}
    maxwell::TM
	F::Vector{Float64}
    MaxwellCF(maxwell::Maxwell{N}) where N = new{N,typeof(maxwell)}(maxwell,deepcopy(maxwell.simulation.F))
end

# Nonlinear Eigenvalue Problem
struct MaxwellNEP{N,TM,TNEP} <: AbstractMaxwellNonlinearEigenproblem{N}
    maxwell::TM
    nep::TNEP
    MaxwellNEP(maxwell::Maxwell{N},NEP) where N = new{N,typeof(maxwell),typeof(NEP)}(maxwell,NEP)
end

# include dimensions
foreach(include,dimensions)

# include interfaces
foreach(include,interfaces)

# add convenience constructors for ElectricField
for ep ∈ (MaxwellLEP,MaxwellCF,MaxwellNEP)
	@eval begin
		# """
		# 	ElectricField(::$($ep),[m=1]) -> ElectricField
		# 	ElectricField(::$($ep),values) -> ElectricField
		#
		# arguments can be provided in any order
		# """
		Common.ElectricField(mep::$(ep)) = ElectricField(mep.sim)
		Common.ElectricField(arg,mep::$(ep)) = ElectricField(mep,arg)
		Common.ElectricField(mep::$(ep),arg) = ElectricField(mep.sim,arg)
	end
end


"""
	maxwelleigen(lep, ω; [ky=0, kz=0, verbose=false, kwargs...]) -> ωs, ψs

Linear eigenvalue solve for frequencies `ωs` closest to `ω`.
`ψ` is an ElectricField.
Keyword `verbose` defaults to `false`.
See docs for ArnoldiMethodWrapper for details of `kwargs`.
------------------
	maxwelleigen(cf, ω; [η=0, verbose, lupack, kwargs...]) -> ηs, us

Linear CF eigenvalues closest to `η`

`us` is an `ElectricField`
Keyword `verbose` defaults to `false`.
Keyword `lupack` defaults to `$DEFAULT_LUPACK` and controls which solver is used in the eigensolver.
See docs for ArnoldiMethodWrapper for details.
------------------
    maxwelleigen(nep, ω; [method=contour_beyn, nk=3, display=false, quad_n=100, kwargs...]) -> ω,ψ
"""
function maxwelleigen(ep::AbstractMaxwellLinearEigenproblem, args...;
				eigensolver::Symbol=DEFAULT_LINEAR_EIGENSOLVER, kwargs...)
	if eigensolver==:ArnoldiMethod
		return maxwelleigen_arnoldimethod(ep, args...; kwargs...)
	elseif eigensolver==:Arpack
		return maxwelleigen_arnoldimethod(ep, args...; kwargs...)
	elseif eigensolver==:IterativeSolvers
		return maxwelleigen_iterativesolvers(ep, args...; kwargs...)
	elseif eigensolver==:KrylovKit
		return maxwelleigen_krylovkit(ep, args...; kwargs...)
	else
		throw("unrecognized eigensolver $eigensolver")
	end
end
function maxwelleigen(ep::AbstractMaxwellNonlinearEigenproblem, args...;
				eigensolver::Symbol=DEFAULT_NONLINEAR_EIGENSOLVER, kwargs...)
	if eigensolver==:NonlinearEigenproblems
		return maxwelleigen_nonlineareigenproblems(ep, args...; kwargs...)
	else
		throw("unrecognized eigensolver $eigensolver")
	end
end


################################################################################
# linear
"""
    MaxwellLEP(::Simulation{1}; ky=0, kz=0) -> lep
"""
MaxwellLEP(args...; kwargs...) = MaxwellLEP(Maxwell(args...; kwargs...))

"""
    lep(ω) -> ∇×∇×, ε+χ(ω)
	lep(ω,ωs,ψs) -> ∇×∇×, ε+χ(ω,ωs,ψs)

`ωs` is an array of frequencies, `ψs` is an ElectricField containing one field for
each element of `ωs`
"""
function (lep::MaxwellLEP)(ω, args...)
    lep.maxwell(ω, args...)
    return lep.maxwell.D², lep.maxwell.αεpχ
end

function Base.getproperty(lep::MaxwellLEP, sym::Symbol)
    if Base.sym_in(sym,(:M,:maxwell))
        return getfield(lep,:maxwell)
    else
        return getproperty(getfield(lep,:maxwell),sym)
    end
end

Base.propertynames(lep::MaxwellLEP, private=false) = private ? fieldnames(MaxwellLEP) : propertynames(lep.maxwell)


################################################################################
# CF
"""
    MaxwellCF(::Simulation{1}; ky=0, kz=0) -> lep
"""
MaxwellCF(args...; kwargs...) = MaxwellCF(Maxwell(args...; kwargs...))

"""
	cf(ω) -> ∇×∇× - εω² , Fω²
	cf(ω,ωs,ψs) -> ∇×∇× - εω² , Fω²/(1 + Σᵢ|γ(ωsᵢ)ψs(i)|²)

`ωs` is an array of frequencies, `ψs` is an ElectricField containing one field for
each element of `ωs`
"""
function (cf::MaxwellCF)(ω,args...)
	cf.maxwell(ω,args...)
	computeF!(cf)
	ω² = ω^2
    rows = rowvals(cf.maxwell.A)
    vals = nonzeros(cf.maxwell.A)
    _, n = size(cf.maxwell.A)
    for i ∈ 1:n
        col = i
        for j ∈ nzrange(cf.maxwell.A, i)
            row = rows[j]
            vals[j] = cf.maxwell.D²[row,col] - ω²*cf.maxwell.αε[row,col]
        end
    end
    return cf.maxwell.A, spdiagm(0=>repeat(ω^2*cf.F,3))
end

# no modification to F if no background field
computeF!(cf::MaxwellCF,ω) = nothing
# modify F (not in simulation, but in cf object), according to background field
function computeF!(cf::MaxwellCF, ω, ωs::Array, ψs::ElectricField)
	χs = flatten(map(d->d.χ,cf.simulation.domains),Common.Dispersions.AbstractDispersion)
    χ = flatten(map(x->susceptability(x,ω),χs),ComplexF64)
	domain_inds = flatten(map(d->ntuple(i->d,length(cf.simulation.domains[d].χ)),ntuple(identity,length(cf.simulation.domains))))
	for j ∈ eachindex(χs)
		if typeof(χs[j])<:TwoLevelSystem
			chi	= χs[j].chi
			for i ∈ eachindex(cf.F)
				if domain_inds[j] == Common.Simulations.which_domain(cf.simulation.domains,cf.simulation.x[i])
					cf.F[i] = cf.simulation.F[i]*real(chi[i]/χ[j])
				end
			end
		end
	end
	return nothing
end

function Base.getproperty(cf::MaxwellCF,sym::Symbol)
    if Base.sym_in(sym,(:M,:maxwell))
        return getfield(cf,:maxwell)
	elseif sym==:F
		return getfield(cf,:F)
    else
        return getproperty(getfield(cf,:maxwell),sym)
    end
end

Base.propertynames(cf::MaxwellCF,private=false) = private ? fieldnames(MaxwellCF) : (propertynames(cf.maxwell)...,:F)


################################################################################
# nonlinear ep

# see dimensional files for dimension-specific constructors

"""
	nep(ω) -> A

full maxwell operator `A` at frequency `ω` for nonlinear eigenproblem `nep`
"""
(mnep::MaxwellNEP)(ω) = compute_Mder(mnep.nep,ω)

function Base.getproperty(nep::MaxwellNEP,sym::Symbol)
    if Base.sym_in(sym,(:M,:maxwell))
        return getfield(nep,:maxwell)
    elseif sym == :nep
        return getfield(nep,:nep)
    else
        return getproperty(getfield(nep,:maxwell),sym)
    end
end

Base.propertynames(nep::MaxwellNEP,private=false) = private ? fieldnames(MaxwellNEP) : (:nep, propertynames(nep.maxwell)...)


################################################################################
# Pretty Printing

function Base.show(io::IO,::MaxwellLEP{N}) where N
    print(io,N,"D ")
    printstyled(io,"MaxwellLEP ",color=PRINTED_COLOR_DARK)
    print(io,"Linear Eigenproblem ")
    printstyled(io,"(pass to ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"maxwelleigen",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", or call w/ args ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ω",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,",[",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ωs",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,",",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ψs",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,"]->",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"A",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,",",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"B",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", where ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"Aψ=ω²Bψ",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
end

function Base.show(io::IO,::MaxwellCF{N}) where N
    print(io,N,"D ")
    printstyled(io,"MaxwellCF ",color=PRINTED_COLOR_DARK)
    print(io,"CF Eigenproblem ")
    printstyled(io,"(pass to ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"maxwelleigen",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", or call w/ args ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ω",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,",[",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ωs",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,",",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ψs",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,"]->",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"A",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,",",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"B",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", where ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"Au=ηBu",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
end

function Base.show(io::IO,::MaxwellNEP{N}) where N
    print(io,N,"D ")
    printstyled(io,"MaxwellNEP ",color=PRINTED_COLOR_DARK)
    print(io,"Nonlinear Eigenproblem for use with ")
    printstyled(io,"(pass to ")
	printstyled(io,"maxwelleigen",color=PRINTED_COLOR_VARIABLE)
	printstyled(io,", or NEP solver like ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"iar_chebyshev",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", or call w/ arg ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ω",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", to evaluate operator)",color=PRINTED_COLOR_INSTRUCTION)
end


end # module

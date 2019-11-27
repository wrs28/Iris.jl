"""
	module Spectral

For the computation of various linear and nonlinear eigenvalue problems.

	*LEP: *L*inear *E*igenvalue *P*roblem, i.e. no bulk dispersion, fixed boundary conditions
		(can use with PMLs/cPMLs for resonance or RSM calculations)
	*CF: *C*onstant *F*lux states, i.e. at a given frequency, find eigen-susceptibilities and
		associated eigenfunctions. Can use either matched boundaries or PMLs. The CF problem
		is a particular kind of linear eigenvalue problem, though here it is treated separately.
	*NEP: *N*onlinear *E*igenvalue *P*roblem, i.e. with bulk dispersion, open boundary conditions,
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
	"Interfaces/KrylovKit.jl",
	"Interfaces/LinearAlgebra.jl",
	"Interfaces/NonlinearEigenproblems.jl"
	)

using ..Common
using SparseArrays
import ..Common.INDEX_OFFSET

abstract type AbstractMaxwellEigenproblem{N} end
abstract type AbstractMaxwellLinearEigenproblem{N} <: AbstractMaxwellEigenproblem{N} end
abstract type AbstractMaxwellNonlinearEigenproblem{N} <: AbstractMaxwellEigenproblem{N} end

# Linear Eigenvalue Problem
struct MaxwellLEP{N,TM} <: AbstractMaxwellLinearEigenproblem{N}
    maxwell::TM
	αεpFχ::SparseMatrixCSC{ComplexF64,Int}
	saturated::Ref{Bool}
    MaxwellLEP(maxwell::Maxwell{N}) where N = new{N,typeof(maxwell)}(maxwell,maxwell.αεpFχ,Ref(false))
end

# Linear Eigenvalue CF Problem
struct MaxwellCF{N,TM} <: AbstractMaxwellLinearEigenproblem{N}
    maxwell::TM
	F::SparseMatrixCSC{Float64,Int}
	saturated::Ref{Bool}
    MaxwellCF(maxwell::Maxwell{N}) where N = new{N,typeof(maxwell)}(maxwell,spdiagm(0=>repeat(maxwell.simulation.F,3)),Ref(false))
end

# Nonlinear Eigenvalue Problem
struct MaxwellNEP{N,TM,TNEP} <: AbstractMaxwellNonlinearEigenproblem{N}
    maxwell::TM
    nep::TNEP
	Fs::Vector{SparseMatrixCSC{ComplexF64,Int}}
	saturated::Ref{Bool}
    MaxwellNEP(maxwell::Maxwell{N},NEP,F_inds) where N = new{N,typeof(maxwell),typeof(NEP)}(maxwell,NEP,NEP.A[F_inds],Ref(false))
end

# include dimensions
foreach(include,dimensions)

# add convenience constructors for ElectricField
for ep ∈ (MaxwellLEP,MaxwellCF,MaxwellNEP)
	@eval begin
		# """
		# 	ElectricField(::$($ep),[m=1]) -> ElectricField
		# 	ElectricField(::$($ep),values) -> ElectricField
		#
		# arguments can be provided in any order
		# """
		Common.ElectricField(mep::$(ep)) = ElectricField(mep.simulation)
		Common.ElectricField(arg,mep::$(ep)) = ElectricField(mep,arg)
		Common.ElectricField(mep::$(ep),arg) = ElectricField(mep.simulation,arg)
	end
end

import ..Common.DEFAULT_LINEAR_EIGENSOLVER
import ..Common.DEFAULT_NONLINEAR_EIGENSOLVER

doc_lep = "
	maxwelleigen(lep, Ω, [ωs, ψs; kwargs...]) -> ω, ψ

Linear eigenvalues `ω` closest to `Ω`, eigenfields `ψ::ElectricField`.

Optional `ωs::Vector` and `ψs::ElectricField` must be provided together, and define
nonlinearity-inducing background fields.

Keywords:"

doc_cf = "
	maxwelleigen(cf, Ω, [ωs, ψs; η=0, kwargs...]) -> H, u

Linear eigen-susceptibilities `H` closest to `η`, eigenfields `u::ElectricField`.
Optional `ωs::Vector` and `ψs::ElectricField` must be provided together, and define
nonlinearity-inducing background fields.

Keywords:
`η` Target CF eigenvalue (`0`);"

doc_nep = "
	maxwelleigen(nep, Ω, [ωs, ψs; kwargs...]) -> ω, ψ

Nonlinear eigenvalues `ω` closest to `Ω`, eigenfields `ψ::ElectricField`.
Optional `ωs::Vector` and `ψs::ElectricField` must be provided together, and define
nonlinearity-inducing background fields.

Keywords:"

# include interfaces
foreach(include,interfaces)

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
	lep.saturated[] = isempty(args) ? false : true
    return lep.maxwell.D², lep.maxwell.αεpFχ
end

function Base.getproperty(lep::MaxwellLEP, sym::Symbol)
	if Base.sym_in(sym,propertynames(getfield(lep,:maxwell)))
        return getproperty(getfield(lep,:maxwell),sym)
	else
		return getfield(lep,sym)
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
	cf.maxwell(ω,args...) # initialize
	computeF!(cf.F,cf.simulation,ω,args...)
	cf.saturated[] = isempty(args) ? false : true
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
	return cf.maxwell.A, ω²*cf.F
end

function Base.getproperty(cf::MaxwellCF,sym::Symbol)
	if Base.sym_in(sym,propertynames(getfield(cf,:maxwell)))
        return getproperty(getfield(cf,:maxwell),sym)
    else
		return getfield(cf,sym)
	end
end

Base.propertynames(cf::MaxwellCF,private=false) = private ? fieldnames(MaxwellCF) : (propertynames(cf.maxwell)...,:F)


################################################################################
# nonlinear ep

# see dimensional files for dimension-specific constructors

"""
	nep(ω) -> A

	nep(ω,ωs,ψs) -> A

full maxwell operator `A` at frequency `ω` for nonlinear eigenproblem `nep`
"""
function (mnep::MaxwellNEP)(ω,args...)
	if !isempty(args)
		mnep.maxwell(ω,args...) #
		mnep.saturated[] = true
	else
		mnep.saturated[] = false
	end
	computeF!(mnep.Fs,mnep.simulation,ω,args...)
	return compute_Mder(mnep.nep,ω)
end

function Base.getproperty(nep::MaxwellNEP,sym::Symbol)
	if Base.sym_in(sym,propertynames(getfield(nep,:maxwell)))
        return getproperty(getfield(nep,:maxwell),sym)
    else
		return getfield(nep,sym)
	end
end

Base.propertynames(nep::MaxwellNEP,private=false) = private ? fieldnames(MaxwellNEP) : (:nep, propertynames(nep.maxwell)...)


################################################################################
# saturate F

# no modification to F if no background field
function computeF!(F::SparseMatrixCSC,sim::Simulation,ω)
	for i ∈ 1:3length(sim)
		j = mod1(i,length(sim))
		if iszero(sim.nondispersive_domain_indices[j])
			F.nzval[i] = 0
		else
			F.nzval[i] = sim.F[j]
		end
	end
	return nothing
end
function computeF!(F::Vector{T},sim::Simulation,ω) where T<:AbstractSparseMatrix
	for d ∈ eachindex(F)
		for i ∈ 1:3length(sim)
			j = mod1(i,length(sim))
			if sim.nondispersive_domain_indices[j] == d
				F[d].nzval[i] = sim.F[j]
			else
				F[d].nzval[i] = 0
			end
		end
	end
	return nothing
end

# modify F (not in simulation, but in cf object), according to background field
function computeF!(F::SparseMatrixCSC, sim::Simulation, ω, ωs::Array, ψs::ElectricField)
	for i ∈ 1:3length(sim)
		j = mod1(i,length(sim))
		if typeof(sim.χ[j])<:TwoLevelSystem
			χ::TwoLevelSystem = sim.χ[j]
			if iszero(χ.D0)
				F[d].nzval[i] = 0
			else
				if iszero(sim.nondispersive_domain_indices[j])
					F[d].nzval[i] = 0
				else
					F[d].nzval[i] = sim.F[j]*χ.hp1⁻¹[i]
				end
			end
		end
	end
	return nothing
end
function computeF!(F::Vector{T}, sim::Simulation, ω, ωs::Array, ψs::ElectricField) where T<:AbstractSparseMatrix
	for d ∈ eachindex(F)
		for i ∈ 1:3length(sim)
			j = mod1(i,length(sim))
			if typeof(sim.χ[j])<:TwoLevelSystem
				χ::TwoLevelSystem = sim.χ[j]
				if iszero(χ.D0)
					F[d].nzval[i] = 0
				else
					if sim.nondispersive_domain_indices[j] == d
						F[d].nzval[i] = sim.F[j]*χ.hp1⁻¹[i]
					else
						F[d].nzval[i] = 0
					end
				end
			end
		end
	end
	return nothing
end

################################################################################

function orthogonalize!(u::ElectricField,sim::Simulation{1},η,B,ky,kz)
	for i ∈ eachindex(η)
		for j ∈ 1:(i-1)
			if abs(2(η[i]-η[j])/(η[i]+η[j])) < 1e-3
				δ = sum(u[:,i].*u[:,j].*diag(B))*sim.dx
				if abs(δ)>.1
					if iszero(ky) && iszero(kz)
						u.z[:,j] .= 0
						u.y[:,i] .= 0
					else
						β = 1/sqrt(1-δ^2)
						α = -δ*β
						u[:,j] = α*u[:,i] + β*u[:,j]
					end
				end
			end
		end
	end
	return nothing
end


################################################################################
# Pretty Printing
import ..Common.PRINTED_COLOR_DARK
import ..Common.PRINTED_COLOR_VARIABLE
import ..Common.PRINTED_COLOR_INSTRUCTION
import ..Common.PRINTED_COLOR_WARN

function Base.show(io::IO,lep::MaxwellLEP{N}) where N
    print(io,N,"D ")
	lep.saturated[] ? printstyled(io,"saturated ",color=PRINTED_COLOR_WARN) : print(io,"unsaturated ")
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

function Base.show(io::IO,cf::MaxwellCF{N}) where N
    print(io,N,"D ")
	cf.saturated[] ? printstyled(io,"saturated ",color=PRINTED_COLOR_WARN) : print(io,"unsaturated ")
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

function Base.show(io::IO,nep::MaxwellNEP{N}) where N
    print(io,N,"D ")
	nep.saturated[] ? printstyled(io,"saturated ",color=PRINTED_COLOR_WARN) : print(io,"unsaturated ")
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

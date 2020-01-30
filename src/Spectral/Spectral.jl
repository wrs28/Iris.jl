"""
Various linear and nonlinear eigenvalue problems.

  * LEP: *L*inear *E*igenvalue *P*roblem, i.e. no bulk dispersion, fixed boundary conditions (can use with PMLs/cPMLs for resonance or RSM calculations)
  * CF: *C*onstant *F*lux states, i.e. at a given frequency, find eigen-susceptibilities and associated eigenfunctions. Can use either matched boundaries or PMLs. The CF problem is a particular kind of linear eigenvalue problem, though here it is treated separately.
  * NEP: *N*onlinear *E*igenvalue *P*roblem, i.e. with bulk dispersion, open boundary conditions, or both. Used to compute resonances/RSMs, no need for PMLs.

API: [`MaxwellLEP`](@ref), [`MaxwellCF`](@ref), [`MaxwellNEP`](@ref), for use with [`maxwelleigen`](@ref) (high-level),
or directly with [`ArnoldiMethod`](https://github.com/haampie/ArnoldiMethod.jl),
[`Arpack`](https://github.com/JuliaLinearAlgebra/Arpack.jl),
[`IterativeSolvers`](https://github.com/JuliaMath/IterativeSolvers.jl),
[`KrylovKit`](https://github.com/Jutho/KrylovKit.jl),
[`LinearAlgebra`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/),
[`NonlinearEigenproblems`](https://github.com/nep-pack/NonlinearEigenproblems.jl) (low-level).
"""
module Spectral

export HelmholtzLEP
export HelmholtzCF
export HelmholtzNEP
export MaxwellLEP
export MaxwellCF
export MaxwellNEP
export helmholtzeigen
export maxwelleigen

using ..Common
using NonlinearEigenproblems

dimensions = (
		"1D/Spectral1D.jl",
        # "2D/Spectral2D.jl",
        # "3D/Spectral3D.jl"
        )

interfaces = (
	# "Interfaces/ArnoldiMethod.jl",
	"Interfaces/Arpack.jl",
	# "Interfaces/KrylovKit.jl",
	"Interfaces/LinearAlgebra.jl",
	"Interfaces/NonlinearEigenproblems.jl"
	)

# defin abstract containers

"""
	AbstractEigenproblem{N}
"""
abstract type AbstractEigenproblem{N} end

"""
	AbstractLinearEigenproblem{N} <: AbstractEigenproblem{N}
"""
abstract type AbstractLinearEigenproblem{N} <: AbstractEigenproblem{N} end

"""
	AbstractCFEigenproblem{N} <: AbstractEigenproblem{N}
"""
abstract type AbstractCFEigenproblem{N} <: AbstractEigenproblem{N} end

"""
	AbstractNonlinearEigenproblem{N} <: AbstractEigenproblem{N}
"""
abstract type AbstractNonlinearEigenproblem{N} <: AbstractEigenproblem{N} end

helmholtzeigen() = nothing
maxwelleigen() = nothing

################################################################################
# Define Eigenproblem Structures for Maxwell & Helmholtz

using SparseArrays

LEPs = (:Helmholtz,:Maxwell)
leps = (:helmholtz,:maxwell)
rep = (1,3) # num of times to repeat F for helmholtz v. maxwell
for (L,l,r) ∈ zip(LEPs,leps,rep)
	@eval begin
		"""
			$($(string(L,"LEP"))){N}

		`N`-dimensional $($(string(L))) Linear Eigenvalue Problem
		"""
		struct $(Symbol(L,"LEP")){N,T} <: AbstractLinearEigenproblem{N}
		    $l::T
			αεpFχ::SparseMatrixCSC{ComplexF64,Int}
			saturated::Ref{Bool}

		    $(Symbol(L,"LEP"))($l::$L{N}) where N = new{N,typeof($l)}($l,$l.αεpFχ,Ref(false))
		end

		"""
			$($(string(L,"CF"))){N}

		`N`-dimensional $($(string(L))) Constant Flux problem
		"""
		struct $(Symbol(L,"CF")){N,T} <: AbstractCFEigenproblem{N}
		    $l::T
			F::SparseMatrixCSC{Float64,Int}
			saturated::Ref{Bool}
		    $(Symbol(L,"CF"))($l::$L{N}) where N = new{N,typeof($l)}($l,spdiagm(0=>repeat($l.simulation.F,$r)),Ref(false))
		end

		"""
			$($(string(L,"NEP"))){N}

		`N`-dimensional $($(string(L))) Constant Flux problem
		"""
		struct $(Symbol(L,"NEP")){N,T,TNEP} <: AbstractNonlinearEigenproblem{N}
		    $l::T
		    nep::TNEP
			Fs::Vector{SparseMatrixCSC{ComplexF64,Int}}
			saturated::Ref{Bool}

		    $(Symbol(L,"NEP"))($l::$L{N},NEP,F_inds) where N = new{N,typeof($l),typeof(NEP)}($l,NEP,NEP.A[F_inds],Ref(false))
		end
	end
end

################################################################################

# include dimensions
foreach(include,dimensions)

################################################################################
# add convenience constructors for VectorFields

# ScalarFields for Helmholtz
for ep ∈ (HelmholtzLEP,HelmholtzCF,HelmholtzNEP)
	@eval begin
		Common.VectorField(hep::$(ep)) = ScalarField(hep)
		Common.VectorField(arg,hep::$(ep)) = ScalarField(hep,arg)
		Common.VectorField(hep::$(ep),arg) = ScalarField(hep,arg)
		Common.ScalarField(hep::$(ep)) = ScalarField(hep.simulation.x)
		Common.ScalarField(arg,hep::$(ep)) = ScalarField(hep.simulation.x,arg)
		Common.ScalarField(hep::$(ep),arg) = ScalarField(hep.simulation.x,arg)
	end
end

# ElectricFields for Maxwell
for ep ∈ (MaxwellLEP,MaxwellCF,MaxwellNEP)
# 	@eval begin
# 		# """
# 		# 	ElectricField(::$($ep),[m=1]) -> ElectricField
# 		# 	ElectricField(::$($ep),values) -> ElectricField
# 		#
# 		# arguments can be provided in any order
# 		# """
# 		Common.ElectricField(mep::$(ep)) = ElectricField(mep.simulation)
# 		Common.ElectricField(arg,mep::$(ep)) = ElectricField(mep.simulation,arg)
# 		Common.ElectricField(mep::$(ep),arg) = ElectricField(mep.simulation,arg)
# 	end
end

################################################################################
# LEP

"""
    HelmholtzLEP(sim) -> lep

`sim` is a `Simulation`, `lep` can be passed to `helmholtzeigen`
"""
HelmholtzLEP(args...; kwargs...) = HelmholtzLEP(Helmholtz(args...; kwargs...))

"""
    MaxwellLEP(sim; ky=0, kz=0) -> lep

`sim` is a `Simulation`, `lep` can be passed to `maxwelleigen`
"""
MaxwellLEP(args...; kwargs...) = MaxwellLEP(Maxwell(args...; kwargs...))

# generate Ax=λBx matrices, with or without saturation for Helmholtz
"""
    lep(ω) -> ∇², -ε-χ(ω)

	lep(ω,ωs,ψs) -> ∇², -ε-χ(ω,ωs,ψs)

`ωs` is an array of frequencies, `ψs` is a `ScalarField` containing one field for
each element of `ωs`
"""
function (lep::HelmholtzLEP)(ω, args...)
    lep.helmholtz(ω, args...)
	lep.saturated[] = isempty(args) ? false : true
    return lep.helmholtz.D², -lep.helmholtz.αεpFχ
end

# generate Ax=λBx matrices, with or without saturation for Maxwell
"""
    lep(ω) -> ∇×∇×, ε+χ(ω)

	lep(ω,ωs,ψs) -> ∇×∇×, ε+χ(ω,ωs,ψs)

`ωs` is an array of frequencies, `ψs` is an `ElectricField` containing one field for
each element of `ωs`
"""
function (lep::MaxwellLEP)(ω, args...)
    lep.maxwell(ω, args...)
	lep.saturated[] = isempty(args) ? false : true
    return lep.maxwell.D², lep.maxwell.αεpFχ
end


################################################################################
# CF

"""
    HelmholtzCF(::Simulation) -> lep
"""
HelmholtzCF(args...; kwargs...) = HelmholtzCF(Helmholtz(args...; kwargs...))
"""
    MaxwellCF(::Simulation; ky=0, kz=0) -> lep
"""
MaxwellCF(args...; kwargs...) = MaxwellCF(Maxwell(args...; kwargs...))

"""
	cf(ω) -> ∇² + εω² , -Fω²

	cf(ω,ωs,ψs) -> ∇² + εω² , -Fω²/(1 + Σᵢ|γ(ωsᵢ)ψs(i)|²)

`ωs` is an array of frequencies, `ψs` is a ScalarField containing one field for
each element of `ωs`
"""
function (cf::HelmholtzCF)(ω,args...)
	cf.helmholtz(ω,args...) # initialize
	_computeF!(cf.F,cf.simulation,ω,args...)
	cf.saturated[] = isempty(args) ? false : true
	ω² = ω^2
    rows = rowvals(cf.helmholtz.A)
    vals = nonzeros(cf.helmholtz.A)
    _, n = size(cf.helmholtz.A)
    for i ∈ 1:n
        col = i
        for j ∈ nzrange(cf.helmholtz.A, i)
            row = rows[j]
            vals[j] = cf.helmholtz.D²[row,col] + ω²*cf.helmholtz.αε[row,col]
        end
    end
	return cf.helmholtz.A, -ω²*cf.F
end

"""
	cf(ω) -> ∇×∇× - εω² , Fω²

	cf(ω,ωs,ψs) -> ∇×∇× - εω² , Fω²/(1 + Σᵢ|γ(ωsᵢ)ψs(i)|²)

`ωs` is an array of frequencies, `ψs` is an ElectricField containing one field for
each element of `ωs`
"""
function (cf::MaxwellCF)(ω,args...)
	cf.maxwell(ω,args...) # initialize
	_computeF!(cf.F,cf.simulation,ω,args...)
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


################################################################################
# NEP
# see dimensional files for dimension-specific constructors

"""
	nep(ω) -> A

	nep(ω,ωs,ψs) -> A

Full helmholtz operator `A` at frequency `ω` for nonlinear eigenproblem `nep`.
Problem saturated by `ωs::Vector` and `ψs::ScalarField`, if provided.
"""
function (mnep::HelmholtzNEP)(ω,args...)
	# saturate if extra args are provided
	if !isempty(args)
		mnep.helmholtz(ω,args...)
		mnep.saturated[] = true
	else
		mnep.saturated[] = false
	end
	# compute F, including saturation
	_computeF!(mnep.Fs, mnep.simulation, ω, args...)
	return compute_Mder(mnep.nep,ω)
end

"""
	nep(ω) -> A

	nep(ω,ωs,ψs) -> A

Full maxwell operator `A` at frequency `ω` for nonlinear eigenproblem `nep`.
Problem saturated by `ωs::Vector` and `ψs::ScalarField`, if provided.
"""
function (mnep::MaxwellNEP)(ω,args...)
	# saturate if extra args are provided
	if !isempty(args)
		mnep.maxwell(ω,args...) #
		mnep.saturated[] = true
	else
		mnep.saturated[] = false
	end
	# compute F, including saturation
	_computeF!(mnep.Fs, mnep.simulation, ω, args...)
	return compute_Mder(mnep.nep,ω)
end

################################################################################
# getproperty and propertynames for LEP, CF, NEP

for (L,l) ∈ zip(LEPs,leps)
	for ep ∈ (:LEP,:CF,:NEP)
		@eval begin
			function Base.getproperty(lep::$(Symbol(L,ep)), sym::Symbol)
				if Base.sym_in(sym,propertynames(getfield(lep,Symbol($(string(l))))))
			        return getproperty(getfield(lep,Symbol($(string(l)))),sym)
				else
					return getfield(lep,sym)
			    end
			end

			Base.propertynames(lep::$(Symbol(L,ep)), private=false) = private ? fieldnames($(Symbol(L,ep))) : propertynames(lep.$l)
		end
	end
end

################################################################################
# saturate F

# CF: No saturation, reset to sim.F
function _computeF!(F::SparseMatrixCSC,sim::Simulation,ω)
	for i ∈ 1:size(F,1)
		j = mod1(i,length(sim))
		# either field is zero if no dispersive domain, or RESET to sim.F value
		if iszero(sim.dispersive_domain_indices[j])
			F.nzval[i] = 0
		else
			F.nzval[i] = sim.F[j]
		end
	end
	return nothing
end

# NEP: No saturation, reset to sim.Fs
function _computeF!(F::Vector{T},sim::Simulation,ω) where T<:AbstractSparseMatrix
	for d ∈ eachindex(F)
		for i ∈ 1:size(F[1],1)
			j = mod1(i,length(sim))
			if sim.dispersive_domain_indices[j] == d
				F[d].nzval[i] = sim.Fs[d][j]
			else
				F[d].nzval[i] = 0
			end
		end
	end
	return nothing
end

# CF: Saturate F (in cf.F, not in sim.F) according to background fields
function _computeF!(F::SparseMatrixCSC, sim::Simulation, ω, ωs::Array, ψs::VectorField)
	for i ∈ 1:size(F,1)
		j = mod1(i,length(sim))
		if typeof(sim.χ[j])<:TwoLevelSystem
			χ::TwoLevelSystem = sim.χ[j]
			if iszero(χ.D0)
				F.nzval[i] = 0
			else
				if iszero(sim.dispersive_domain_indices[j])
					F.nzval[i] = 0
				else
					F.nzval[i] = sim.F[j]*χ.hp1⁻¹[i]
				end
			end
		end
	end
	return nothing
end

# NEP: Saturate F (not in sim.F) according to background fields
function _computeF!(F::Vector{T}, sim::Simulation, ω, ωs::Array, ψs::VectorField) where T<:AbstractSparseMatrix
	for d ∈ eachindex(F)
		for i ∈ 1:size(F[1],1)
			j = mod1(i,length(sim))
			if typeof(sim.χ[j])<:TwoLevelSystem
				χ::TwoLevelSystem = sim.χ[j]
				if iszero(χ.D0)
					F[d].nzval[i] = 0
				else
					if sim.dispersive_domain_indices[j] == d
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
# ORTHOGONALIZE

using LinearAlgebra
import ..Common.ORTHOGONALIZE_OVERLAP_THRESHOLD

"""
	orthogonalize!(ψ,sim,η,B,ky,kz,[ind])

Check for nearly degenerate eigenpairs, if overlap integral exceeds $ORTHOGONALIZE_OVERLAP_THRESHOLD,
orthoganlize and renormalize.
"""
function orthogonalize!(u::ScalarField{N},sim::Simulation{N},η,B,args...) where N
	for i ∈ eachindex(η)
		for j ∈ 1:(i-1)
			if abs(2(η[i]-η[j])/(η[i]+η[j])) < 1e-3
				δ = sum(u[:,i].*u[:,j].*diag(B))*_measure(sim)
				if abs(δ) > ORTHOGONALIZE_OVERLAP_THRESHOLD
					β = 1/sqrt(1-δ^2)
					α = -δ*β
					u[:,j] = α*u[:,i] + β*u[:,j]
				end
			end
		end
	end
	return nothing
end

function orthogonalize!(u::ElectricField{N},sim::Simulation{N},η,B,ky,kz,args...) where N
	for i ∈ eachindex(η)
		for j ∈ 1:(i-1)
			if abs(2(η[i]-η[j])/(η[i]+η[j])) < 1e-3
				δ = sum(u[:,i].*u[:,j].*diag(B))*_measure(sim)
				if abs(δ) > ORTHOGONALIZE_OVERLAP_THRESHOLD
					if all(iszero,(ky,kz))
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
	if all(iszero,(ky,kz)) normalize!(sim,u,B,args...) end
	return nothing
end

################################################################################
# Pretty Printing

import ..Common.PRINTED_COLOR_DARK
import ..Common.PRINTED_COLOR_VARIABLE
import ..Common.PRINTED_COLOR_INSTRUCTION
import ..Common.PRINTED_COLOR_WARN

# LEP
function Base.show(io::IO,lep::AbstractLinearEigenproblem{N}) where N
    print(io,N,"D ")
	lep.saturated[] ? printstyled(io,"saturated ",color=PRINTED_COLOR_WARN) : nothing#print(io,"unsaturated ")
	if typeof(lep) <: MaxwellLEP
	    printstyled(io,"MaxwellLEP ",color=PRINTED_COLOR_DARK)
	elseif typeof(lep) <: HelmholtzLEP
	    printstyled(io,"HelmholtzLEP ",color=PRINTED_COLOR_DARK)
	end
    printstyled(io,"(pass to ",color=PRINTED_COLOR_INSTRUCTION)
	if typeof(lep) <: MaxwellLEP
	    printstyled(io,"maxwelleigen",color=PRINTED_COLOR_VARIABLE)
	elseif typeof(lep) <: HelmholtzLEP
	    printstyled(io,"helmholtzeigen",color=PRINTED_COLOR_VARIABLE)
	end
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

# CF
function Base.show(io::IO,cf::AbstractCFEigenproblem{N}) where N
    print(io,N,"D ")
	cf.saturated[] ? printstyled(io,"saturated ",color=PRINTED_COLOR_WARN) : nothing#print(io,"unsaturated ")
	if typeof(cf) <: MaxwellCF
		printstyled(io,"MaxwellCF ",color=PRINTED_COLOR_DARK)
	elseif typeof(cf) <: HelmholtzCF
		printstyled(io,"HelmholtzCF ",color=PRINTED_COLOR_DARK)
	end
    printstyled(io,"(pass to ",color=PRINTED_COLOR_INSTRUCTION)
	if typeof(cf) <: MaxwellCF
		printstyled(io,"maxwelleigen",color=PRINTED_COLOR_VARIABLE)
	elseif typeof(cf) <: HelmholtzCF
		printstyled(io,"helmholtzeigen",color=PRINTED_COLOR_VARIABLE)
	end
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

# NEP
function Base.show(io::IO,nep::AbstractNonlinearEigenproblem{N}) where N
    print(io,N,"D ")
	nep.saturated[] ? printstyled(io,"saturated ",color=PRINTED_COLOR_WARN) : nothing#print(io,"unsaturated ")
	if typeof(nep) <: MaxwellNEP
	    printstyled(io,"MaxwellNEP ",color=PRINTED_COLOR_DARK)
	elseif typeof(nep) <: HelmholtzNEP
		printstyled(io,"HelmholtzNEP ",color=PRINTED_COLOR_DARK)
	end
    printstyled(io,"(pass to ", color=PRINTED_COLOR_INSTRUCTION)
	if typeof(nep) <: MaxwellNEP
	    printstyled(io,"maxwelleigen",color=PRINTED_COLOR_VARIABLE)
	elseif typeof(nep) <: HelmholtzNEP
		printstyled(io,"helmholtzeigen",color=PRINTED_COLOR_VARIABLE)
	end
	printstyled(io,", or NEP solver like ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"iar_chebyshev",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", or call w/ arg ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ω",color=PRINTED_COLOR_VARIABLE)
	printstyled(io,",[",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ωs",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,",",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ψs",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,"]->",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"A",color=PRINTED_COLOR_VARIABLE)
	printstyled(io,", where ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"Aψ=0",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
end


################################################################################
# DOCSTRINGS
# to be paired with default interface for `helmholtzeigen` and `maxwelleigen`

DOC_LEP_H = "
	helmholtzeigen(lep, Ω, [ωs, ψs; kwargs...]) -> ω, ψ

Linear eigenvalues `ω` closest to `Ω`, eigenfields `ψ::ScalarField`.

Optional `ωs::Vector` and `ψs::ScalarField` must be provided together, and define
nonlinearity-inducing background fields.

Keywords:"


DOC_LEP_M = "
	maxwelleigen(lep, Ω, [ωs, ψs; kwargs...]) -> ω, ψ

Linear eigenvalues `ω` closest to `Ω`, eigenfields `ψ::ElectricField`.

Optional `ωs::Vector` and `ψs::ElectricField` must be provided together, and define
nonlinearity-inducing background fields.

Keywords:"


DOC_CF_H = "
	helmholtzeigen(cf, Ω, [ωs, ψs; η=0, kwargs...]) -> H, u

Linear eigen-susceptibilities `H` closest to `η`, eigenfields `u::ScalarField`.
Optional `ωs::Vector` and `ψs::ScalarField` must be provided together, and define
nonlinearity-inducing background fields.

Keywords:
`η` Target CF eigenvalue (`0`);"


DOC_CF_M = "
	maxwelleigen(cf, Ω, [ωs, ψs; η=0, kwargs...]) -> H, u

Linear eigen-susceptibilities `H` closest to `η`, eigenfields `u::ElectricField`.
Optional `ωs::Vector` and `ψs::ElectricField` must be provided together, and define
nonlinearity-inducing background fields.

Keywords:
`η` Target CF eigenvalue (`0`);"


DOC_NEP_H = "
	helmholtzeigen(nep, Ω, [ωs, ψs; kwargs...]) -> ω, ψ

Nonlinear eigenvalues `ω` closest to `Ω`, eigenfields `ψ::ScalarField`.
Optional `ωs::Vector` and `ψs::ScalarField` must be provided together, and define
nonlinearity-inducing background fields.

Keywords:"


DOC_NEP_M = "
	maxwelleigen(nep, Ω, [ωs, ψs; kwargs...]) -> ω, ψ

Nonlinear eigenvalues `ω` closest to `Ω`, eigenfields `ψ::ElectricField`.
Optional `ωs::Vector` and `ψs::ElectricField` must be provided together, and define
nonlinearity-inducing background fields.

Keywords:"

################################################################################
# Interfaces

import ..Common.DEFAULT_LINEAR_EIGENSOLVER
import ..Common.DEFAULT_NONLINEAR_EIGENSOLVER
import ..Common.INDEX_OFFSET

foreach(include,interfaces)

end # module

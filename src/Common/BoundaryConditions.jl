"""
for constructing boundary conditions used in the definition of a Boundary
"""
module BoundaryConditions

export AbstractBC
export AbstractComplexBC
export AbstractRealBC
export getside
export noBC
export DirichletBC
export NeumannBC
export FloquetBC
export MatchedBC
export LocalBC
export NonLocalBC
export BCLocality
export BCHermiticity
export HermitianBC
export NonHermitianBC

import ..PRINTED_COLOR_DARK
import ..getside

"""
	AbstractBC{SIDE}
"""
abstract type AbstractBC{SIDE} end
"""
	AbstractRealBC{SIDE} <: AbstractBC{SIDE}
"""
abstract type AbstractRealBC{SIDE} <: AbstractBC{SIDE} end
"""
	AbstractComplexBC{SIDE} <: AbstractBC{SIDE}
"""
abstract type AbstractComplexBC{SIDE} <: AbstractBC{SIDE} end

"""
	getside(::AbstractBC) -> side
"""
getside(::AbstractBC{SIDE}) where SIDE = SIDE

"""
	noBC{SIDE} <: AbstractRealBC{SIDE}
"""
struct noBC{SIDE} <: AbstractRealBC{SIDE} end

"""
	DirichletBC{SIDE} <: AbstractRealBC{SIDE}
"""
struct DirichletBC{SIDE} <: AbstractRealBC{SIDE} end

"""
	NeumannBC{SIDE} <: AbstractRealBC{SIDE}
"""
struct NeumannBC{SIDE} <: AbstractRealBC{SIDE} end

"""
	FloquetBC{SIDE} <: AbstractComplexBC{SIDE}

NOT YET IMPLEMENTED
"""
struct FloquetBC{SIDE} <: AbstractComplexBC{SIDE} end

"""
	MatchedBC{SIDE} <: AbstractComplexBC{SIDE}

fields: `in`, `out`
"""
struct MatchedBC{SIDE} <: AbstractComplexBC{SIDE}
	in::Array{Int,1}
	out::Array{Int,1}
end
"""
	MatchedBC{SIDE}([in, out]) -> mbc

`in::Vector` are input channels specified by integer
`out::Vector` are output channels specified by integer
"""
MatchedBC{SIDE}(;in::Array{Int,1}=Int[],out::Array{Int,1}=Int[]) where SIDE = MatchedBC{SIDE}(in,out)

"""
	LocalBC{SIDE}
"""
struct LocalBC{SIDE} end

"""
	LocalBC{SIDE}
"""
struct NonLocalBC{SIDE} end

BCLocality(::Union{noBC{SIDE},DirichletBC{SIDE},NeumannBC{SIDE}}) where SIDE = LocalBC{SIDE}
"""
	BCLocality(bc) -> locality

get locality trait of a boundary condition
"""
BCLocality(::AbstractBC{SIDE}) where SIDE = NonLocalBC{SIDE}

"""
	HermitianBC{SIDE}
"""
struct HermitianBC{SIDE} end

"""
	NonHermitianBC{SIDE}
"""
struct NonHermitianBC{SIDE} end

BCHermiticity(::Union{noBC{SIDE},DirichletBC{SIDE},NeumannBC{SIDE},FloquetBC{SIDE}}) where SIDE = HermitianBC{SIDE}
"""
	BCHermiticity(bc) -> hermiticity

get hermiticity trait of a boundary condition
"""
BCHermiticity(::AbstractBC{SIDE}) where SIDE = NonHermitianBC{SIDE}

for TBC ∈ (noBC,DirichletBC,NeumannBC,FloquetBC)
	@eval begin
		function Base.show(io::IO,::$TBC{SIDE}) where SIDE
			printstyled(io,$TBC,color=PRINTED_COLOR_DARK)
			!get(io,:compact,false) ? print(io," SIDE=",SIDE) : nothing
		end
	end
end

# Pretty Printing
function Base.show(io::IO,mbc::MatchedBC{SIDE}) where SIDE
	printstyled(io,"MatchedBC",color=PRINTED_COLOR_DARK)
	!get(io,:compact,false) ? print(io," SIDE=",SIDE,";") : nothing
	print(IOContext(io,:typeinfo=>Array{Int})," in:", mbc.in)
	print(io,"/out:", mbc.out)
end

Base.conj(bc::AbstractRealBC) = bc

Base.conj(bc::MatchedBC{SIDE}) where SIDE = MatchedBC{SIDE}(in=bc.out,out=bc.in)

end # module

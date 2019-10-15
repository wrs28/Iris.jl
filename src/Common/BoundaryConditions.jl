"""
	module BoundaryConditions
"""
module BoundaryConditions

export noBC
export DirichletBC
export NeumannBC
export FloquetBC
export MatchedBC
export AbstractBC
export AbstractRealBC
export AbstractComplexBC
export LocalBC
export NonLocalBC
export BCLocality

abstract type AbstractBC{SIDE} end
abstract type AbstractRealBC{SIDE} <: AbstractBC{SIDE} end
abstract type AbstractComplexBC{SIDE} <: AbstractBC{SIDE} end

struct noBC{SIDE} <: AbstractRealBC{SIDE} end
struct DirichletBC{SIDE} <: AbstractRealBC{SIDE} end
struct NeumannBC{SIDE} <: AbstractRealBC{SIDE} end
struct FloquetBC{SIDE} <: AbstractComplexBC{SIDE} end
struct MatchedBC{SIDE} <: AbstractComplexBC{SIDE}
	in::Array{Int,1}
	out::Array{Int,1}
end
MatchedBC{SIDE}(;in::Array{Int,1}=Int[],out::Array{Int,1}=Int[]) where SIDE = MatchedBC{SIDE}(in,out)

struct LocalBC{SIDE} end
struct NonLocalBC{SIDE} end

BCLocality(::Union{noBC{SIDE},DirichletBC{SIDE},NeumannBC{SIDE}}) where SIDE = LocalBC{SIDE}
BCLocality(::AbstractBC{SIDE}) where SIDE = NonLocalBC{SIDE}

for TBC ∈ (noBC,DirichletBC,NeumannBC,FloquetBC)
	@eval Base.show(io::IO,::$TBC{SIDE}) where SIDE = print(io,$TBC," (SIDE=",SIDE,")")
end
function Base.show(io::IO,mbc::MatchedBC{SIDE}) where SIDE
	print(io,"MatchedBC (SIDE=",SIDE,";")
	print(io," in:", mbc.in)
	print(io," out:", mbc.out)
	print(io,")")
end

Base.conj(bc::AbstractRealBC) = bc
Base.conj(bc::MatchedBC{SIDE}) where SIDE = MatchedBC{SIDE}(in=bc.out,out=bc.in)

"""
	struct noBC{SIDE}
------
	noBC{SIDE}(...) -> bc
"""
noBC

"""
	struct DirichletBC{SIDE}
------
	DirichletBC{SIDE}(depth) -> bc
"""
DirichletBC

"""
	struct NeumannBC{SIDE}
------
	NeumannBC{SIDE}(depth) -> bc
"""
NeumannBC

"""
	struct FloquetBC{SIDE}
------
	FloquetBC{SIDE}(depth) -> bc
"""
FloquetBC

"""
	struct MatchedBC{SIDE}
------
	MatchedBC{SIDE}(depth) -> bc
"""
MatchedBC

end # module

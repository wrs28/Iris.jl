"""
	module BoundaryConditions
"""
module BoundaryConditions

export get_side
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

import ..PRINTED_COLOR_DARK
import ..get_side

abstract type AbstractBC{SIDE} end
abstract type AbstractRealBC{SIDE} <: AbstractBC{SIDE} end
abstract type AbstractComplexBC{SIDE} <: AbstractBC{SIDE} end

get_side(::AbstractBC{SIDE}) where SIDE = SIDE

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
	@eval begin
		function Base.show(io::IO,::$TBC{SIDE}) where SIDE
			printstyled(io,$TBC,color=PRINTED_COLOR_DARK)
			!get(io,:compact,false) ? print(io," SIDE=",SIDE) : nothing
		end
	end
end
function Base.show(io::IO,mbc::MatchedBC{SIDE}) where SIDE
	printstyled(io,"MatchedBC",color=PRINTED_COLOR_DARK)
	!get(io,:compact,false) ? print(io," SIDE=",SIDE,";") : nothing
	print(IOContext(io,:typeinfo=>Array{Int})," in:", mbc.in)
	print(io,"/out:", mbc.out)
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

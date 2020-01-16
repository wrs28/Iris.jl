"""
for constructing boundary layers and evaluating [`PML`](@ref)/[`cPML`](@ref) conductivity
"""
module BoundaryLayers

export PML
export cPML
export noBL
export AbstractBL
export AbstractRealBL
export AbstractComplexBL
export getside

# PML params used in Boundaries
import ..EXTINCTION
import ..SCALING_ANGLE
import ..BL_DEPTH

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK

using ..Points
using Formatting

"""
	AbstractBL{SIDE}
"""
abstract type AbstractBL{SIDE} end
"""
	AbstractRealBL{SIDE} <: AbstractBL{SIDE}
"""
abstract type AbstractRealBL{SIDE} <: AbstractBL{SIDE} end
"""
	AbstractComplexBL{SIDE} <: AbstractBL{SIDE}
"""
abstract type AbstractComplexBL{SIDE} <: AbstractBL{SIDE} end

"""
	getside(::AbstractBL) = side
"""
getside(::AbstractBL{SIDE}) where SIDE = SIDE

for (TBL,ATBL) ∈ zip((:PML,:cPML,:noBL),
					 (:AbstractComplexBL,:AbstractComplexBL,:AbstractRealBL))
	@eval begin
		"""
			$($TBL){SIDE}

		Fields: `depth`, `start`, `stop`
		"""
		mutable struct $TBL{SIDE} <: $ATBL{SIDE}
			depth::Float64
		    start::Float64
		    stop::Float64

			function $TBL{SIDE}(depth,start,stop) where SIDE
				if !($TBL <: noBL)
					depth>0 || throw(ArgumentError("cannot have $($TBL) with 0 depth"))
				end
				return new{SIDE}(depth,start,stop)
			end
		end

		if $TBL == :noBL
			noBL{SIDE}(depth::Real=BL_DEPTH) where SIDE = noBL{SIDE}(0,NaN,NaN)
		else
			"""
				$($TBL){SIDE}(depth) -> bl

			$($TBL) along dimension `SIDE` with `depth`. The location of the layer determined
			later by associating it with a shape in a call to `Boundary`.

			see also: [`Boundary`](@ref), [`PML`](@ref), [`cPML`](@ref), [`noBL`](@ref)
			"""
			$TBL{SIDE}(depth::Real=BL_DEPTH) where SIDE = $TBL{SIDE}(depth,NaN,NaN)
		end

		function Base.show(io::IO,tbl::$TBL{SIDE}) where SIDE
			printstyled(io,$TBL,color=PRINTED_COLOR_DARK)
			!get(io,:compact,false) ? print(io," SIDE=",SIDE) : nothing
			if $TBL!==noBL
				print(io," depth=")
				printstyled(io,fmt("1.1f",tbl.depth),color=PRINTED_COLOR_NUMBER)
			end
		end

		(tbl::$TBL)(x::Real,y...) = tbl(Point(x,y...))

		(bl::$TBL{SIDE})(;depth=bl.depth, start=bl.start, stop=bl.stop) where SIDE = $TBL{SIDE}(depth,start,stop)
	end
end
# noBL{SIDE}(depth::Real=BL_DEPTH) where SIDE = noBL{SIDE}(0,NaN,NaN)

"""
	(::PML)(p::Point) -> σ

evaluate the PML conductivity `σ` at the point `p`.
"""
(pml::PML)(p::Point) = _conductivity_profile(p[1+(getside(pml)-1)÷2],pml.start,pml.stop)

"""
	(::cPML)(p::Point) -> σ

evaluate the cPML conductivity `σ` at the point `p`.
"""
(cpml::cPML)(p::Point) = _conj_conductivity_profile(p[1+(getside(cpml)-1)÷2],cpml.start,cpml.stop)

"""
	(::noBL)(p::Point) -> 0.0
"""
(nbl::noBL)(p::Point) = 0.0

"""
	conj(pml) -> cpml
"""
Base.conj(pml::PML{SIDE}) where SIDE = cPML{SIDE}(pml.depth,pml.start,pml.stop)
"""
	conj(cpml) -> pml
"""
Base.conj(cpml::cPML{SIDE}) where SIDE = PML{SIDE}(cpml.depth,cpml.start,cpml.stop)
"""
	conj(nobl) -> nobl
"""
Base.conj(nbl::noBL) = nbl

# used inside PML. This is the place to modify the PML profile
"""
	_conductivity_profile(z,start,stop) -> σ

PML conductivity `σ` at point `z`
"""
function _conductivity_profile(z::Real,start::Real,stop::Real)
	if (start≤stop && start≤z) || (start≥stop && start≥z)
		depth = abs(stop-start)
	    β = 2*sqrt(depth)
	    α = -(1/4)*exp(complex(0,SCALING_ANGLE))*log(EXTINCTION)/depth/((2exp(β)-β)/(2β*expm1(β))-1/β^2)
        u = abs(z-start)/depth
        return α*u*expm1(β*u)/expm1(β)
    else
        return complex(0.0)
    end
end

"""
	_integrated_conductivity_profile(z,start,stop) -> σ

integrated PML conductivity `σ` at point `z`
"""
function _integrated_conductivity_profile(z::Real,start::Real,stop::Real)
	if (start≤stop && start≤z) || (start≥stop && start≥z)
		depth = abs(stop-start)
		β = 2*sqrt(depth)
		α = -(1/4)*exp(complex(0,SCALING_ANGLE))*log(EXTINCTION)/depth/((2exp(β)-β)/(2β*expm1(β))-1/β^2)
		u = abs(z-start)/depth
		return (α*depth/expm1(β))*(u*exp(β*u)-(1+1/β)*expm1(β*u)/β)
	else
		return complex(0.0)
	end
end


# used inside cPML
"""
	_conj_conductivity_profile(z,start,stop) -> σ

conjugate PML conductivity `σ` at point `z`
"""
_conj_conductivity_profile(z,start,stop) = -conj(_conductivity_profile(z,start,stop))

"""
	_conj_integrated_conductivity_profile(z,start,stop) -> σ

integrated conjugate PML conductivity `σ` at point `z`
"""
_conj_integrated_conductivity_profile(z,start,stop) = -conj(_integrated_conductivity_profile(z,start,stop))

end # module

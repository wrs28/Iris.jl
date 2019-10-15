"""
	module BoundaryLayers
"""
module BoundaryLayers

export PML
export cPML
export noBL
export AbstractBL

using ...Defaults
using ..Points
using Formatting

abstract type AbstractBL{SIDE} end
abstract type AbstractRealBL{SIDE} <: AbstractBL{SIDE} end
abstract type AbstractComplexBL{SIDE} <: AbstractBL{SIDE} end

for (TBL,ATBL) ∈ zip((:PML,:cPML,:noBL),
					 (:AbstractComplexBL,:AbstractComplexBL,:AbstractRealBL))
	@eval begin
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

		$TBL{SIDE}(depth) where SIDE = $TBL{SIDE}(depth,NaN,NaN)

		function Base.show(io::IO,tbl::$TBL{SIDE}) where SIDE
			print(io,$TBL," (SIDE=",SIDE, ", depth=",fmt("1.1f",tbl.depth),")")
		end

		(tbl::$TBL)(x::Real,y...) = tbl(Point(x,y...))

		(bl::$TBL{SIDE})(;depth=bl.depth, start=bl.start, stop=bl.stop) where SIDE = $TBL{SIDE}(depth,start,stop)
	end
end
noBL{SIDE}() where SIDE = noBL{SIDE}(0,NaN,NaN)

(bl::PML{SIDE})(p::Point) where {SIDE} = conductivity_profile(p[ceil(Int,SIDE/2)],bl.start,bl.stop)
(bl::cPML{SIDE})(p::Point) where {SIDE} = conj_conductivity_profile(p[ceil(Int,SIDE/2)],bl.start,bl.stop)
(bl::noBL)(p::Point) = 0.0


"""
	struct PML{SIDE}(depth) -> pml

absorbing PML along dimension `SIDE` with `depth`. The location of the PML is later determined
by associating it with a shape in a call to `Boundary`.

see also: [`Boundary`](@ref), [`cPML`](@ref), [`noBL`](@ref)

---

	(::PML)(x,y...) -> σ

evaluate the PML conductivity `σ` at the point `x,y...` in local frame.
"""
PML

"""
	struct cPML{SIDE}(depth) -> pml

amplifying conjugate PML along dimension `SIDE` with `depth`. The location of the cPML is later determined
by associating it with a shape in a call to `Boundary`.

see also: [`Boundary`](@ref), [`PML`](@ref), [`noBL`](@ref)

---

	(::cPML)(x,y...) -> σ

evaluate the PML conductivity `σ` at the point `x,y...` in local frame
"""
cPML

"""
	struct noBL{SIDE}(depth) -> nobl

lack of boundary layer along dimension `SIDE`. `depth` is not used, but it there for consistency.
The location is later determined by associating it with a shape in a call to `Boundary`.

see also: [`Boundary`](@ref), [`PML`](@ref), [`cPML`](@ref)

---

	(::noBL)(x,y...) -> σ

returns `σ` = 0. This is here for consistency.
"""
noBL


Base.conj(pml::PML{SIDE}) where SIDE = cPML{SIDE}(pml.depth,pml.start,pml.stop)
Base.conj(cpml::cPML{SIDE}) where SIDE = PML{SIDE}(cpml.depth,cpml.start,cpml.stop)
Base.conj(nbl::noBL) = nbl
"""
	conj(pml) -> cpml
	conj(cpml) -> pml
	conj(nobl) -> nobl

conjugate a boundary layer
"""
conj


# used inside PML. This is the place to modify the PML profile
function conductivity_profile(z::Real,start::Real,stop::Real)
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
function integrated_conductivity_profile(z::Real,start::Real,stop::Real)
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
conj_conductivity_profile(z,start,stop) = -conj(conductivity_profile(z,start,stop))
conj_integrated_conductivity_profile(z,start,stop) = -conj(integrated_conductivity_profile(z,start,stop))


end # module

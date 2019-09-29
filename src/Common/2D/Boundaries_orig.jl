"""
	module Boundaries
"""
module Boundaries

export Boundary,
noBC,
DirichletBC,
NeumannBC,
FloquetBC,
MatchedBC,
AbstractBC,
AbstractLocalBC,
AbstractRealBC,
AbstractComplexBC,
PML,
cPML,
noBL,
AbstractBL,
AbstractRealBL,
AbstractRealBC

using ...Defaults
using ..Shapes
using RecipesBase

abstract type AbstractBC end
abstract type AbstractRealBC <: AbstractBC end
abstract type AbstractComplexBC <: AbstractBC end
abstract type AbstractLocalBC <: AbstractRealBC end

struct noBC <: AbstractBC end
struct DirichletBC <: AbstractLocalBC end
struct NeumannBC <: AbstractLocalBC end
struct FloquetBC <: AbstractComplexBC end
struct MatchedBC <: AbstractComplexBC end

# because displaying struct with no fields is broken
Base.show(io::IO,bc::noBC) = print(io,"noBC")
Base.show(io::IO,bc::DirichletBC) = print(io,"Dirichlet")
Base.show(io::IO,bc::NeumannBC) = print(io,"Neumann")
Base.show(io::IO,bc::FloquetBC) = print(io,"Floquet")
Base.show(io::IO,bc::MatchedBC) = print(io,"Matched")

abstract type AbstractBL end
abstract type AbstractRealBL <: AbstractBL end
abstract type AbstractComplexBL <: AbstractBL end

"""
	struct PML(depth) -> pml

absorbing PML with depth `depth`. The location of the PML is later determined
by associating it with a shape in a call to `Boundary`.

see also: [`Boundary`](@ref), [`cPML`](@ref), [`noBL`](@ref)

---

	(::PML)(x,y) -> σ

evaluate the PML conductivity `σ` at the point `x,y` in local frame.
"""
struct PML <: AbstractComplexBL
	depth::Float64 # depth of PML layer
    startx::Float64
    stopx::Float64
	starty::Float64
    stopy::Float64

	function PML(depth,startx,stopx,starty,stopy)
		@assert depth>0 "cannot have PML with 0 depth"
		return new(depth,startx,stopx,starty,stopy)
	end
    PML(depth::Real) = PML(depth,NaN,NaN,NaN,NaN)
    Base.show(io::IO,pml::PML) = print(io,"PML")
end
"""
	struct cPML(depth) -> pml

amplifying conjugate PML with depth `depth`. The location of the cPML is later determined
by associating it with a shape in a call to `Boundary`.

see also: [`Boundary`](@ref), [`PML`](@ref), [`noBL`](@ref)

---

	(::cPML)(x,y) -> σ

evaluate the PML conductivity `σ` at the point `x,y` in local frame
"""
struct cPML <: AbstractComplexBL
	depth::Float64 # depth of PML layer
    startx::Float64
    stopx::Float64
	starty::Float64
    stopy::Float64

    function cPML(depth,startx,stopx,starty,stopy)
        @assert depth>0 "cannot have cPML with 0 depth"
        return new(depth,startx,stopx,starty,stopy)
    end
	cPML(depth::Real) = cPML(depth,NaN,NaN,NaN,NaN)
    Base.show(io::IO,cpml::cPML) = print(io,"cPML")
end
"""
	struct noBL(depth) -> nobl

lack of boundary layer. `depth` is not used, but it there for consistency.
The location is later determined by associating it with a shape in a call to `Boundary`.

see also: [`Boundary`](@ref), [`PML`](@ref), [`cPML`](@ref)

---

	(::noBL)(x,y) -> σ

returns `σ` = 0. This is here for consistency.
"""
struct noBL <: AbstractRealBL
	depth::Float64 # depth of PML layer
    startx::Float64
    stopx::Float64
	starty::Float64
    stopy::Float64

    noBL(args...) = new(0,NaN,NaN,NaN,NaN)
    Base.show(io::IO,nobl::noBL) = print(io,"noBL")
end

TBL = (PML,cPML,noBL)
@eval for tbl ∈ TBL
	function (bl::tbl)(;depth=bl.depth, startx=bl.startx, stopx=bl.stopx, starty=bl.starty, stopy=bl.stopy)
		return tbl(depth,startx,stopx,starty,stopy)
	end
	function (bl::tbl)(x::Number,y::Number)
		if tbl!==noBL
			@assert  isfinite(bl.stopx-bl.startx) || isfinite(bl.stopy-bl.starty) "boundary layer not well defined. check startx and stopx fields (and y)"
			z = isfinite(bl.stopx-bl.startx) ? x : y
			start = isfinite(bl.stopx-bl.startx) ? bl.startx : bl.starty
			stop = isfinite(bl.stopx-bl.startx) ? bl.stopx : bl.stopy
			if tbl<:PML
				return conductivity_profile(z,start,stop)
			elseif tbl<:cPML
				return conj_conductivity_profile(z,start,stop)
			end
		else
			return complex(0.0)
		end
	end
end


"""
	conj(::PML) -> cpml::cPML
	conj(::cPML) -> pml::PML
	conj(::noBL) -> nobl::noBL

conjugate a boundary layer
"""
Base.conj(pml::PML) = cPML(pml.depth,pml.startx,pml.stopx,pml.starty,pml.stopy)
Base.conj(cpml::cPML) = PML(cpml.depth,cpml.startx,cpml.stopx,cpml.starty,cpml.stopy)
Base.conj(nbl::noBL) = noBL()


# used inside PML. This is the place to modify the PML profile
function conductivity_profile(z::Real,start::Real,stop::Real)
	if (start≤stop && start≤z) || (start≥stop && start≥z)
		depth = abs(stop-start)
	    β = 2*sqrt(depth)
	    α = -(1/4)*exp(complex(0,SCALING_ANGLE))*log(EXTINCTION)/depth/((2exp(β)-β)/(2β*expm1(β))-1/β^2)
        u = abs(z-start)/depth
        σ = α*u*expm1(β*u)/expm1(β)
    else
        σ = complex(0.0)
    end
    return σ
end

# used inside cPML
conj_conductivity_profile(x::Real,start::Real,stop::Real) = -conj(conductivity_profile(x,start,stop))

function integrated_conductivity_profile(z::Real,start::Real,stop::Real)
	if (start≤stop && start≤z) || (start≥stop && start≥z)
		depth = abs(stop-start)
		β = 2*sqrt(depth)
		α = -(1/4)*exp(complex(0,SCALING_ANGLE))*log(EXTINCTION)/depth/((2exp(β)-β)/(2β*expm1(β))-1/β^2)
		u = abs(z-start)/depth
		Σ = (α*depth/expm1(β))*(u*exp(β*u)-(1+1/β)*expm1(β*u)/β)
	else
		Σ = complex(0.0)
	end
	return Σ
end

conj_integrated_conductivity_profile(x::Real,start::Real,stop::Real) = -conj(integrated_conductivity_profile(x,start,stop))


"""
	struct Boundary

	Boundary(shape::Shape,bcs,bls) -> bnd

`bcs` is either tuple of boundary conditions, or it is a single boundary condition which applies to all sides.
`bls` is either tuple of boundary layers, or it is a single boundary layer which applies to all sides.

if fewer boundary conditions or layers are supplied than the number of sides of
`shape`, it assumes `noBC` and `noBL`
"""
struct Boundary{TS,TBC,TBL}
    shape::TS
    bcs::TBC
    bls::TBL

	function Boundary(shape::AbstractShape{N},bcs::AbstractBC=noBC(),bls::AbstractBL=noBL()) where N
		bcs = fill(bcs,N)
		bls = fill(bls,N)
		return Boundary(shape,(bcs...,),(bls...,))
	end
	function Boundary(shape::AbstractShape{N},bcs::AbstractBC,bls) where N
		bcs = fill(bcs,N)
		return Boundary(shape,(bcs...,),bls)
	end
	function Boundary(shape::AbstractShape{N},bcs,bls::AbstractBL=noBL()) where N
		bls = fill(bls,N)
		return Boundary(shape,bcs,(bls...,))
	end
    function Boundary(shape::AbstractShape{N},
        bcs::Tuple,
        bls::Tuple
        ) where {N}

		if length(bcs)<N
			@warn "fewer boundary conditions provided ($(length(bcs))) than shape requires ($N). rest assumed to be noBC"
		elseif length(bcs)>N
			@warn "more boundary conditions provided ($(length(bcs))) than shape requires ($N). ignoring all bc's after $N"
			bcs = bcs[1:N]
		end
		if length(bls)>N
			@warn "more boundary layers provided ($(length(bls))) than shape requires ($N). ignoring all bl's after $N"
			bls = bls[1:N]
		end

		bcs_type = promote_type(typeof.(bcs)...)
		length(bcs)≤N ? bcs_type=promote_type(bcs_type,noBC) : nothing
		BCS = bcs_type[bcs...]
		for i ∈ 1:(N-length(bcs))
			append!(BCS,[noBC()])
		end

        bls_type = promote_type(typeof.(bls)...)
        length(bls)≤N ? bls_type=promote_type(bls_type,noBL) : nothing
        BLS = bls_type[bls...]
		for i ∈ 1:(N-length(bls))
	        append!(BLS,[noBL()])
        end

		fix_bls_by_shape!(BLS,shape)

		BLS = (BLS...,)
        BCS = (BCS...,)

        return new{typeof(shape),typeof(BCS),typeof(BLS)}(shape,BCS,BLS)
    end

	function Base.show(io::IO,bnd::Boundary)
		N = get_number_of_sides(bnd.shape)
		print(io,"Boundary: \n")
		print(io,"\tshape: $(bnd.shape)")
		for i ∈ 1:N
			println(io)
			print(io,"\tside ",i,": ")
			print(io,bnd.bcs[i],"/",bnd.bls[i])
		end
	end
end


"""
	conj(::Boundary) -> bnd

Make a new `Boundary` object with conjugated boundary layers.
"""
Base.conj(bnd::Boundary) = Boundary(bnd.shape,bnd.bcs,conj.(bnd.bls))


# set location of boundary layer. Valid only for squares, rectangles, circle,
# since for now these are the only shapes supporting PMLs/cPMLs.
fix_bls_by_shape!(bls,s::Square) = fix_bls_by_shape!(bls,Rectangle(s))
function fix_bls_by_shape!(bls,shape::Rectangle)
	stop = shape.x0-shape.a/2
	start = stop + bls[1].depth
	bls[1] = bls[1](startx=start, stopx=stop)

	stop = shape.x0+shape.a/2
	start = stop - bls[2].depth
	bls[2] = bls[2](startx=start, stopx=stop)

	stop = shape.y0-shape.b/2
	start = stop + bls[3].depth
	bls[3] = bls[3](starty=start, stopy=stop)

	stop = shape.y0+shape.b/2
	start = stop - bls[4].depth
	bls[4] = bls[4](starty=start, stopy=stop)
end
function fix_bls_by_shape!(bls,c::Circle)
	stop=c.R
	start=stop-bls[1].depth
	bls[1] = bls[1](startx=start,stopx=stop)
end
fix_bls_by_shape!(bls,shape::AbstractShape) = nothing


"""
	plot(::Boundary)

PML's plotted in blue patches (b/c absorbing, so cold), cPML's in red (b/c emitting, so hot)

DirichletBC are solid black lines, NeumannBC dashed, FloquetBC dashdot, MatchedBC/noBC nothing
"""
@recipe function f(bnd::Boundary)
	aspect_ratio --> 1
	legend --> false

	@series begin
		fillcolor --> BOUNDARY_COLOR
		color --> BOUNDARY_COLOR
		markershape := :none
		bnd.shape
	end
	@series (bnd,1)
	@series (bnd,1,1)
end
# plot boundary layers
@recipe function f(bnd::Boundary,_1::Int)
	shape = bnd.shape
	cosθ = shape.cosθ
	sinθ = shape.sinθ
	typeof(shape)<:Square ? shape = Rectangle(shape) : nothing
	for i ∈ eachindex(bnd.bls)
		bls = bnd.bls[i]
		depth = bnd.bls[i].depth
		if typeof(shape)<:Rectangle
			if i==1
				x = -(shape.a-depth)/2
				y = 0
				bl_shape = Rectangle(depth,shape.b,shape.x0+cosθ*x-sinθ*y,shape.y0+sinθ*x+cosθ*y,shape.θ)
			elseif i==4
				x = 0
				y = (shape.b-depth)/2
				bl_shape = Rectangle(shape.a,depth,shape.x0+cosθ*x-sinθ*y,shape.y0+sinθ*x+cosθ*y,shape.θ)
			elseif i==2
				x = +(shape.a-depth)/2
				y = 0
				bl_shape = Rectangle(depth,shape.b,shape.x0+cosθ*x-sinθ*y,shape.y0+sinθ*x+cosθ*y,shape.θ)
			else
				x = 0
				y = -(shape.b-depth)/2
				bl_shape = Rectangle(shape.a,depth,shape.x0+cosθ*x-sinθ*y,shape.y0+sinθ*x+cosθ*y,shape.θ)
			end
		elseif typeof(shape)<:Circle
			bls = bnd.bls[1]
			depth = bls.depth
			bl_shape = Annulus(shape.R-depth,shape.R,shape.x0,shape.y0)
		end
		if typeof(bls)<:AbstractComplexBL
			@series begin
				alpha --> 0
				typeof(bls)<:PML ? fillcolor --> BOUNDARY_PML_COLOR : fillcolor --> BOUNDARY_CPML_COLOR
				typeof(bls)<:PML ? color --> BOUNDARY_PML_COLOR : color --> BOUNDARY_CPML_COLOR
				bl_shape
			end
		end
	end
end
# plot boundary conditions
@recipe function f(bnd::Boundary,_1::Int,_2::Int)
	shape = bnd.shape
	for i ∈ eachindex(bnd.bcs)
		bcs = bnd.bcs[i]
		if typeof(shape)<:Rectangle
			if i==1
				x = [-1,-1]*shape.a/2
				y = [-1,+1]*shape.b/2
			elseif i==2
				x = [-1,+1]*shape.a/2
				y = [+1,+1]*shape.b/2
			elseif i==3
				x = [+1,+1]*shape.a/2
				y = [-1,+1]*shape.b/2
			else
				x = [-1,+1]*shape.a/2
				y = [-1,-1]*shape.b/2
			end
			if typeof(bcs)<:AbstractLocalBC
				@series begin
					if typeof(bcs)<:DirichletBC
						ls --> BOUNDARY_DIRICHLET_LINETYPE
					elseif typeof(bcs)<:NeumannBC
						ls --> BOUNDARY_NEUMANN_LINETYPE
					else
						ls --> BOUNDARY_FLOQUET_MATCHED_LINETYPE
					end
					seriestype --> :line
					markershape := :none
					color --> BOUNDARY_BC_COLOR
					shape.x0 .+ shape.cosθ*x-shape.sinθ*y, shape.y0 .+ shape.sinθ*x+shape.cosθ*y
				end
			end
		else
			if typeof(bcs)<:Union{AbstractLocalBC,FloquetBC}
				@series begin
					if typeof(bcs)<:DirichletBC
						ls --> BOUNDARY_DIRICHLET_LINETYPE
					elseif typeof(bcs)<:NeumannBC
						ls --> BOUNDARY_NEUMANN_LINETYPE
					else
						ls --> BOUNDARY_FLOQUET_MATCHED_LINETYPE
					end
					color --> BOUNDARY_BC_COLOR
					markershape := :none
					fillalpha --> 0
					alpha --> 1
					shape
				end
			end
		end
	end
end


end # module

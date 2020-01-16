# I think this file is done 01/03/2020 WRS

"""
N-dimensional Point containers
"""
module Points

export Point
export Cartesian
export Polar
export Spherical

using Formatting
using LinearAlgebra
using RecipesBase
using StaticArrays

abstract type AbstractCoordinateType end
struct Cartesian <: AbstractCoordinateType end
struct Polar <: AbstractCoordinateType end
struct Spherical <: AbstractCoordinateType end

"""
	Point{N,C}

`N`-Dimensional point of coordinate type `C`
Can be accessed with fields `:x`, `:y`, `:z`, `:s`, `:r`, `:ϕ`, `:θ` and their variants.
"""
struct Point{N,C}
	vector::SVector{N,Float64}

	function Point{N,C}(vec::SVector{N,T}) where {N,C<:AbstractCoordinateType,T}
		!any(isnan.(vec)) || throw(ArgumentError("argument of `Point` $(sum(isnan.(vec))) contains NaNs"))
		all(isfinite.(vec)) || throw(ArgumentError("argument of `Point` contains $(sum(isfinite.(vec))) Infs"))
		new{N,C}(vec)
	end

	Point{Cartesian}(vec::SVector{N,T}) where {N,T<:Real} = Point{N,Cartesian}(vec)

	function Point{Polar}(vec::SVector{N,T}) where {N,T<:Real}
		1≤N≤3 || throw(ArgumentError("`Point{Polar}` requires 1≤N≤3, but N=$N"))
		0≤vec[1] || throw(ArgumentError("`Point{Polar}` requires 0≤r, but r=$(vec[1])"))
		return Point{N,Polar}(vec)
	end
	function Point{Spherical}(vec::SVector{N,T}) where {N,T<:Real}
		1≤N≤3 || throw(ArgumentError("`Point{Spherical}` requires 1≤N≤3, but N=$N"))
		0≤vec[1] || throw(ArgumentError("Point{Spherical}` requires 0≤r, but r=$(vec[1])"))
		if N==3 0≤vec[3]≤π || throw(ArgumentError("`Point{Spherical}` requires 0≤θ≤π, but θ=$(vec[3])")) end
		return Point{N,Spherical}(vec)
	end
end

# default Point constructor makes N-d Cartesian point
"""
	Point(x₁,x₂...,xⱼ) -> point::Point{j,Cartesian}
"""
Point(arg,args...) = Point{Cartesian}(arg,args...)

# Point{C} constructor takes an ntuple or a list of numbers
"""
	Point{C}(x₁,x₂...,xⱼ) -> point::Point{j,C}

`C` is an `AbstractCoordinateType`, one of `Cartesian`, `Polar`, `Spherical`
"""
Point{C}(vec::Vararg{Number,N}) where {N,C} = Point{C}(float.(vec))
Point{C}(vec::NTuple{N}) where {N,C} = Point{C}(SVector{N,Float64}(vec))

# make a new point from `p` with data in `vec`
Point(p::Point) = p
Point(p::Point{N,C},vec::SVector{N}) where {N,C} = Point{C}(vec)
Point{C}(p::Point{N,C},vec::SVector{N}) where {N,C} = Point{C}(vec)

# conversion to Cartesian point
function Cartesian(p::Point{N,C}) where {N,C}
	if C<:Cartesian
		return p
	elseif C<:Polar
		if N==1
			return Point{Cartesian}(p.x)
		elseif N==2
			return Point{Cartesian}(p.x,p.y)
		elseif N==3
			return Point{Cartesian}(p.x,p.y,p.z)
		else
			throw(ArgumentError("conversion from `Polar` to `Cartesian` in dimension N=$N not implemented"))
		end
	elseif C<:Spherical
		if N==1
			return Point{Cartesian}(p.x)
		elseif N==2
			return Point{Cartesian}(p.x,p.z)
		elseif N==3
			return Point{Cartesian}(p.x,p.y,p.z)
		else
			throw(ArgumentError("conversion from `Spherical` to `Cartesian` in dimension N=$N not implemented"))
		end
	else
		throw(ArgumentError("cannot convert `Point{$C}` to `Point{Cartesian}`"))
	end
end
# conversion to Polar point
function Polar(p::Point{N,C}) where {N,C}
	if C<:Polar
		return p
	elseif C<:Cartesian
		if N==1
			return Point{Polar}(p.s)
		elseif N==2
			return Point{Polar}(p.s,p.ϕ)
		elseif N==3
			return Point{Polar}(p.s,p.ϕ,p.z)
		else
			throw(ArgumentError("conversion from `Cartesian` to `Polar` in dimension N=$N not implemented"))
		end
	elseif C<:Spherical
		if N==1
			return Point{Polar}(p.s)
		elseif N==2
			throw(ArgumentError("conversion from `Spherical` to `Polar` in dimension N=2 not implemented"))
		elseif N==3
			return Point{Polar}(p.s,p.ϕ,p.z)
		else
			throw(ArgumentError("conversion from `Spherical` to `Polar` in dimension N=$N not implemented"))
		end
	else
		throw(ArgumentError("cannot convert `Point{$C}` to `Point{Polar}`"))
	end
end
# conversion to Spherical point
function Spherical(p::Point{N,C}) where {N,C}
	if C<:Spherical
		return p
	elseif C<:Cartesian
		if N==1
			return Point{Spherical}(p.r)
		elseif N==2
			return Point{Spherical}(p.r,π/2-p.ϕ)
		elseif N==3
			return Point{Spherical}(p.r,p.θ,p.ϕ)
		else
			throw(ArgumentError("conversion from `Cartesian` to `Spherical` in dimension N=$N not implemented"))
		end
	elseif C<:Polar
		if N==1
			return Point{Spherical}(p.r)
		elseif N==2
			return Point{Spherical}(p.r,p.θ)
		elseif N==3
			return Point{Spherical}(p.r,p.θ,p.ϕ)
		else
			throw(ArgumentError("conversion from `Polar` to `Spherical` in dimension N=$N not implemented"))
		end
	else
		throw(ArgumentError("cannot convert `Point{$C}` to `Point{Polar}`"))
	end
end

# define action of `convert` by using `Cartesian`, `Polar`, `Spherical`
Base.convert(::Type{Point{N,C}}, p::Point) where {N,C} = C(p)
# 1-dim conversion of point to single float
Base.convert(::Type{Float64}, p::Point{1}) = convert(Float64, p.x)

# for use in tests (see defs of norm later)
function Base.isapprox(x::Point, y::Point; atol::Real=0, rtol::Real=sqrt(eps(Float64)), nans::Bool=false)
    x == y || norm(x-y) ≤ max(atol, rtol*max(norm(x), norm(y)))
end

Base.:*(p::Point,a::Number) = a*p
Base.:*(a::Number,p::Point{N,Cartesian}) where N = Point{Cartesian}(a*p.vector)
Base.:*(a::Number,p::Point{1,Polar}) = Point{Polar}(a*p.vec[1])
Base.:*(a::Number,p::Point{2,Polar}) = Point{Polar}(a*p.vec[1],p.vec[2])
Base.:*(a::Number,p::Point{3,Polar}) = Point{Polar}(a*p.vec[1],p.vec[2],a*p.vec[3])
Base.:*(a::Number,p::Point{1,Spherical}) = Point{Spherical}(a*p.vec[1])
Base.:*(a::Number,p::Point{2,Spherical}) = Point{Spherical}(a*p.vec[1],p.vec[2])
Base.:*(a::Number,p::Point{3,Spherical}) = Point{Spherical}(a*p.vec[1],p.vec[2],p.vec[3])

Base.:\(a::Number,p::Point{N}) where N = p/a
Base.:/(p::Point,a::Number) where N = p*(1/a)

Base.:+(p1::Point{N,Cartesian},p2::Point{N,Cartesian}) where N = Point{Cartesian}(p1.vector+p2.vector)
Base.:-(p1::Point{N,Cartesian},p2::Point{N,Cartesian}) where N = Point{Cartesian}(p1.vector-p2.vector)
# differing types defaults to Cartesian
Base.:+(p1::Point{N,C1},p2::Point{N,C2}) where {N,C1,C2} = Cartesian(p1) + Cartesian(p2)
Base.:-(p1::Point{N,C1},p2::Point{N,C2}) where {N,C1,C2} = Cartesian(p1) - Cartesian(p2)
# same types preserves coordinate type
Base.:+(p1::Point{N,C},p2::Point{N,C}) where {N,C} = C(Cartesian(p1) + Cartesian(p2))
Base.:-(p1::Point{N,C},p2::Point{N,C}) where {N,C} = C(Cartesian(p1) - Cartesian(p2))

Base.:-(p::Point) = -1*p
Base.:-(p::Point{1,Polar}) = p
Base.:-(p::Point{2,Polar}) = Point{Polar}(p.s,mod2pi(p.ϕ+π))
Base.:-(p::Point{3,Polar}) = Point{Polar}(p.s,mod2pi(p.ϕ+π),-p.z)
Base.:-(p::Point{1,Spherical}) = p
Base.:-(p::Point{2,Spherical}) = Point{Spherical}(p.r, mod(p.θ+π,π))
Base.:-(p::Point{3,Spherical}) = Point{Spherical}(p.r, mod(p.ϕ+π,π), mod2pi(p.ϕ+π))

Base.:+(p1::Point{1},p2::Number) = Point{Cartesian}(p1.vec .+ p2)
Base.:+(p2::Number, p1::Point{1}) = p1 + p2
Base.:-(p1::Point{1},p2::Number) = +(p1,-p2)
Base.:-(p2::Number,p1::Point{1}) = +(p2,-p1)

LinearAlgebra.norm(p::Point{N}) where N = p.r

Base.ndims(::Point{N}) where N = N
Base.getindex(p::Point,index::Integer) = getfield(p,:vector)[index]
fnames = (:firstindex, :lastindex, :length, :size, :eltype, :iterate)
for fname ∈ fnames @eval Base.$fname(p::Point) = $fname(p.vec) end
Base.iterate(p::Point,state) = iterate(p.vec,state)
Base.IndexStyle(::Point) = IndexStyle(p.vec)

Base.isless(p1::Point{1},p2::Point{1}) = isless(p1.x,p2.x)
Base.isless(p::Point{1},x::Real) = isless(p.x,x)
Base.isless(x::Real,p::Point{1}) = isless(x,p.x)

################################################################################
# Properties

# Cartesian
function Base.getproperty(p::Point{N,Cartesian},sym::Symbol) where N
	if sym == :vec
		return getfield(p,:vector)
	elseif Base.sym_in(sym,(:x,:X))
		return p[1]
	elseif Base.sym_in(sym,(:y,:Y))
		if N≤1
			return 0.0
		else
			return p[2]
		end
	elseif Base.sym_in(sym,(:z,:Z))
		if N≤2
			return 0.0
		else
			return p[3]
		end
	elseif Base.sym_in(sym,(:s,))
		if N==1
			return abs(p[1])
		elseif N≥2
			return hypot(p[1],p[2])
		end
	elseif Base.sym_in(sym,(:r,))
		return norm(getfield(p,:vector))
	elseif Base.sym_in(sym,(:ϕ,:φ,:𝝋,:phi))
		if N==1
			return p[1]>0 ? 0.0 : π
		elseif 2≤N≤3
			return atan(p[2],p[1])
		else
			throw(ErrorException("$sym only defined for 1≤N≤3, but N=$N"))
		end
	elseif Base.sym_in(sym,(:θ,:theta,:ϑ))
		if N≤2
			return π/2
		elseif N==3
			θ = acos(p[3]/norm(getfield(p,:vector)))
			return isnan(θ) ? 0.0 : θ
		else
			throw(ErrorException("$sym only defined for 1≤N≤3, but N=$N"))
		end
	else
		return getfield(p,sym)
	end
end

# Polar
function Base.getproperty(p::Point{N,Polar},sym::Symbol) where N
	if sym == :vec
		return getfield(p,:vector)
	elseif Base.sym_in(sym,(:x,:X))
		if N==1
			return p[1]
		else
			return p[1]*cos(p[2])
		end
	elseif Base.sym_in(sym,(:y,:Y))
		if N==1
			return 0.0
		else
			return p[1]*sin(p[2])
		end
	elseif Base.sym_in(sym,(:z,:Z))
		if N≤2
			return 0.0
		else
			return p[3]
		end
	elseif Base.sym_in(sym,(:s,))
		return p[1]
	elseif Base.sym_in(sym,(:r,))
		if N≤2
			return p[1]
		else
			return hypot(p[1],p[3])
		end
	elseif Base.sym_in(sym,(:ϕ,:φ,:𝝋,:phi))
		if N==1
			return 0.0
		else
			return p[2]
		end
	elseif Base.sym_in(sym,(:θ,:theta,:ϑ))
		if N≤2
			return π/2
		else
			return acos(p[3]/hypot(p[1],p[3]))
		end
	else
		return getfield(p,sym)
	end
end

# Spherical
function Base.getproperty(p::Point{N,Spherical},sym::Symbol) where N
	if sym == :vec
		return getfield(p,:vector)
	elseif Base.sym_in(sym,(:x,:X))
		if N==1
			return p[1]
		elseif N==2
			return p[1]*sin(p[2])
		elseif N==3
			return p[1]*sin(p[2])*cos(p[3])
		else
			throw(ErrorException("N=$N and field :x not consistent with Point{Spherical}"))
		end
	elseif Base.sym_in(sym,(:y,:Y))
		if N≤2
			return 0.0
		else
			return p[1]*sin(p[2])*sin(p[3])
		end
	elseif Base.sym_in(sym,(:z,:Z))
		if N==1
			return 0.0
		else
			return p[1]*cos(p[2])
		end
	elseif Base.sym_in(sym,(:s,))
		if N==1
			return p[1]
		else
			return p[1]*sin(p[2])
		end
	elseif Base.sym_in(sym,(:r,))
		return p[1]
	elseif Base.sym_in(sym,(:θ,:ϑ,:𝜃,:theta))
		if N==1
			return 0.0
		else
			return p[2]
		end
	elseif Base.sym_in(sym,(:ϕ,:φ,:𝝋,:phi))
		if N≤2
			return 0.0
		else
			return p[3]
		end
	else
		return getfield(p,sym)
	end
end

#1D
function Base.propertynames(::Point{1,C},private=false) where C
	if private
		return fieldnames(Point)
	elseif C<:Cartesian
		return (:x, :vector)
	elseif C<:Polar
		return (:s, :vector)
	elseif C<:Spherical
		return (:r, :vector)
	end
end

#2D
function Base.propertynames(::Point{2,C},private=false) where C
	if private
		return fieldnames(Point)
	elseif C<:Cartesian
		return (:x, :y, :vector, :s, :ϕ)
	elseif C<:Polar
		return (:s, :ϕ, :vector, :x, :y)
	elseif C<:Spherical
		return (:r, :θ, :vector, :x, :z)
	end
end

#3D
function Base.propertynames(::Point{3,C},private=false) where C
	if private
		return fieldnames(Point)
	elseif C<:Cartesian
		return (:x, :y, :z, :vector, :r, :ϕ, :θ, :s)
	elseif C<:Polar
		return (:s, :ϕ, :z, :vector, :x, :y, :r, :θ)
	elseif C<:Spherical
		return (:r, :ϕ, :θ, :vector, :x, :y, :z, :s)
	end
end

################################################################################
# Pretty Printing

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK

# 1d (Cartesian)
function Base.show(io::IO,p::Point{1,Cartesian})
	print(io,"1D Cartesian ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": ")
	printstyled(io, fmt("3.3f",p.x),color=PRINTED_COLOR_NUMBER)
end
# 1d (Polar)
function Base.show(io::IO,p::Point{1,Polar})
	print(io,"1D Polar ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": s=")
	printstyled(io, fmt("3.3f",p.s),color=PRINTED_COLOR_NUMBER)
end
# 1d (Spherical)
function Base.show(io::IO,p::Point{1,Spherical})
	print(io,"1D Spherical ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": r=")
	printstyled(io, fmt("3.3f",p.r),color=PRINTED_COLOR_NUMBER)
end

# 2d (Cartesian)
function Base.show(io::IO,p::Point{2,Cartesian})
	print(io,"2D Cartesian ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": x=")
	printstyled(io, fmt("3.3f",p.x),color=PRINTED_COLOR_NUMBER)
	print(io, ", y=")
	printstyled(io, fmt("3.3f",p.y),color=PRINTED_COLOR_NUMBER)
end
# 2d (Polar)
function Base.show(io::IO,p::Point{2,Polar})
	print(io,"2D Polar ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": s=")
	printstyled(io, fmt("3.3f",p.s),color=PRINTED_COLOR_NUMBER)
	print(io, ", ϕ=")
	printstyled(io, fmt("3.3f",180*p.ϕ/π),"°",color=PRINTED_COLOR_NUMBER)
end
# 2d (Spherical)
function Base.show(io::IO,p::Point{2,Spherical})
	print(io,"2D Spherical ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": r=")
	printstyled(io, fmt("3.3f",p.r),color=PRINTED_COLOR_NUMBER)
	print(io, ", θ=")
	printstyled(io, fmt("3.3f",180*p.θ/π),"°",color=PRINTED_COLOR_NUMBER)
end

# 3d (Cartesian)
function Base.show(io::IO,p::Point{3,Cartesian})
	print(io,"3D Cartesian ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": x=")
	printstyled(io, fmt("3.3f",p.x),color=PRINTED_COLOR_NUMBER)
	print(io, ", y=")
	printstyled(io, fmt("3.3f",p.y),color=PRINTED_COLOR_NUMBER)
	print(io, ", z=")
	printstyled(io, fmt("3.3f",p.z),color=PRINTED_COLOR_NUMBER)
end
# 3d (Polar)
function Base.show(io::IO,p::Point{3,Polar})
	print(io,"3D Polar ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": s=")
	printstyled(io, fmt("3.3f",p.s),color=PRINTED_COLOR_NUMBER)
	print(io, ", ϕ=")
	printstyled(io, fmt("3.3f",180*p.ϕ/π),"°",color=PRINTED_COLOR_NUMBER)
	print(io, ", z=")
	printstyled(io, fmt("3.3f",p.z),color=PRINTED_COLOR_NUMBER)
end
# 3d (Spherical)
function Base.show(io::IO,p::Point{3,Spherical})
	print(io,"3D Spherical ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": r=")
	printstyled(io, fmt("3.3f",p.r),color=PRINTED_COLOR_NUMBER)
	print(io, ", θ=")
	printstyled(io, fmt("3.3f",180*p.θ/π),"°",color=PRINTED_COLOR_NUMBER)
	print(io, ", ϕ=")
	printstyled(io, fmt("3.3f",180*p.ϕ/π),"°",color=PRINTED_COLOR_NUMBER)
end

function Base.show(io::IO,p::Point{N,C}) where {N,C}
	print(io, "$(N)D $C ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
end

################################################################################
# Plotting

@recipe function f(p::Point{1})
	seriestype --> :scatter
	legend --> false
	([p.x],[0])
end

@recipe function f(p::Point{2})
	seriestype --> :scatter
	legend --> false
	([p.x],[p.y])
end

@recipe function f(p::Point{3})
	seriestype --> :scatter
	legend --> false
	([p.x],[p.y],[p.z])
end

@recipe function f(p::Array{T}) where T<:Point{1}
	seriestype --> :scatter
	legend --> false
	x = map(p -> p.x, p)
	y = map(p -> 0, p)
	(x,y)
end

@recipe function f(p::Array{T}) where T<:Point{2}
	seriestype --> :scatter
	legend --> false
	x = map(p -> p.x, p)
	y = map(p -> p.y, p)
	(x,y)
end

@recipe function f(p::Array{T}) where T<:Point{3}
	seriestype --> :scatter
	legend --> false
	x = map(p -> p.x, p)
	y = map(p -> p.y, p)
	z = map(p -> p.z, p)
	(x,y,z)
end

end # module

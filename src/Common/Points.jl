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

	function Point{C}(vec::SVector{N,T}) where {C<:AbstractCoordinateType,N,T<:Real}
		C <: AbstractCoordinateType || throw("incorrect coordinate type $C")
		if N==1 C<:Cartesian || throw("N=1 requires coordinate type `Cartesian`") end
		if C<:Polar 2≤N≤3 || throw("type `Polar` requires N>1, but N=$N") end
		if C<:Spherical N==3 || throw("type `Cartesian` requires N>2, but N=$N") end
		return new{N,C}(vec)
	end
end

"""
	Point(x,y...) -> point::Point
"""
Point(x::Real,y...) = Point{Cartesian}(x,y...)

"""
	Point{C}(x,y...) -> point::Point{N,C}

`C` is an `AbstractCoordinateType`, one of `Cartesian`, `Polar`, `Spherical`
"""
Point{C}(x::Real,y...) where C = Point{C}((float(x),float.(y)...))

Point(vec::NTuple{N}) where N = Point{Cartesian}(vec)
Point{C}(vec::NTuple{N}) where {N,C} = Point{C}(SVector{N,Float64}(vec))

Point(p::Point) = p
Point(p::Point{N,C},vec::SVector{N}) where {N,C} = Point{C}(p,vec)
Point{C}(p::Point{N,C},vec::SVector{N}) where {N,C} = Point{C}(p.vec)

function Cartesian(p::Point{N,C}) where {N,C}
	if C<:Cartesian
		return p
	else
		if N==2
			return Point{Cartesian}((p.x,p.y))
		else
			return Point{Cartesian}((p.x,p.y,p.z))
		end
	end
end
function Polar(p::Point{N,C}) where {N,C}
	if C<:Polar
		return p
	else
		if N==2
			return Point{Polar}((p.s,p.ϕ))
		else
			return Point{Polar}((p.s,p.ϕ,p.z))
		end
	end
end
function Spherical(p::Point{N,C}) where {N,C}
	if C<:Spherical
		return p
	else
		return Point{Spherical}((p.r,p.ϕ,p.θ))
	end
end

Base.convert(::Type{Point{N,C}}, p::Point) where {N,C} = C(p)

Base.:*(p::Point,a::Number) = *(a,p)
Base.:*(a::Number,p::Point{N,Cartesian}) where N = Point{Cartesian}(a*p.vector)
Base.:*(a::Number,p::Point{N,Polar}) where N = Point{Polar}((a*p.vec[1],Base.tail(p.vector.data)...))
Base.:*(a::Number,p::Point{N,Spherical}) where N = Point{Spherical}((a*p.vec[1],Base.tail(p.vector.data)...))

Base.:\(a::Number,p::Point{N}) where N = p/a
Base.:/(p::Point{N,Cartesian},a::Number) where N = Point{Cartesian}(p.vector/a)
Base.:/(p::Point{N,Polar},a::Number) where N = Point{Polar}((p.vec[1]/a,Base.tail(p.vector.data)...))
Base.:/(p::Point{N,Spherical},a::Number) where N = Point{Spherical}((p.vec[1]/a,Base.tail(p.vector.data)...))

Base.:+(p1::Point{N,Cartesian},p2::Point{N,Cartesian}) where N = Point{Cartesian}(p1.vec+p2.vector)
# differing types defaults to Cartesian
Base.:+(p1::Point{N,C1},p2::Point{N,C2}) where {N,C1,C2} = Cartesian(p1) + Cartesian(p2)
# same types preserves coordinate type
Base.:+(p1::Point{N,C},p2::Point{N,C}) where {N,C} = C(Cartesian(p1) + Cartesian(p2))

Base.:-(p1::Point,p2::Point) = +(p1,*(-1,p2))

Base.:+(p1::Point{1},p2::Number) = Point{Cartesian}(p1.vec .+ p2)
Base.:+(p2::Number,p1::Point{1}) = +(p1,p2)
Base.:-(p1::Point{1},p2::Number) = +(p1,-p2)
Base.:-(p2::Number,p1::Point{1}) = +(p2,-p1)

LinearAlgebra.norm(p::Point{N,Cartesian}) where N = norm(p.vec)
LinearAlgebra.norm(p::Point{2,Polar}) = p.s
LinearAlgebra.norm(p::Point{3,Polar}) = hypot(p.s,p.z)
LinearAlgebra.norm(p::Point{3,Spherical}) = p.r

Base.convert(::Type{Float64}, p::Point{1}) = convert(Float64, p.x)

Base.ndims(::Point{N}) where N = N
Base.getindex(p::Point,index::Integer) = getfield(p,:vector)[index]
fnames = (:firstindex, :lastindex, :length, :size, :eltype, :iterate)
for fname ∈ fnames @eval Base.$fname(p::Point) = $fname(p.vec) end
Base.iterate(p::Point,state) = iterate(p.vec,state)
Base.IndexStyle(::Point) = IndexStyle(p.vec)

Base.isless(p1::Point{1},p2::Point{1}) = isless(p1.x,p2.x)
Base.isless(p::Point{1},x::Real) = isless(p.x,x)
Base.isless(x::Real,p::Point{1}) = isless(x,p.x)

function Base.getproperty(p::Point{N,Cartesian},sym::Symbol) where N
	if sym == :vec
		return getfield(p,:vector)
	elseif Base.sym_in(sym,(:x,:X))
		N≥1 || throw(ErrorException("N=$N not compatible with $sym"))
		return p[1]
	elseif Base.sym_in(sym,(:y,:Y))
		N≥2 || throw(ErrorException("N=$N not compatible with $sym"))
		return p[2]
	elseif Base.sym_in(sym,(:z,:Z))
		N≥3 || throw(ErrorException("N=$N not compatible with $sym"))
		return p[3]
	elseif Base.sym_in(sym,(:s,))
		N≥2 || throw(ErrorException("N=$N not compatible with $sym"))
		return hypot(p[1],p[2])
	elseif Base.sym_in(sym,(:r,))
		return norm(getfield(p,:vector))
	elseif Base.sym_in(sym,(:ϕ,:φ,:𝝋,:phi))
		2≤N≤3 || throw(ErrorException("$sym only defined for 2≤N≤3, but N=$N"))
		return atan(p[2],p[1])
	elseif Base.sym_in(sym,(:θ,:theta,:ϑ))
		N==3 || throw(ErrorException("$sym only defined for N=3, but N=$N"))
		return acos(p[3]/norm(getfield(p,:vector)))
	else
		return getfield(p,sym)
	end
end
function Base.getproperty(p::Point{N,Polar},sym::Symbol) where N
	if sym == :vec
		return getfield(p,:vector)
	elseif Base.sym_in(sym,(:x,:X))
		N≥1 || throw(ErrorException("N=$N not compatible with $sym"))
		return p[1]*cos(p[2])
	elseif Base.sym_in(sym,(:y,:Y))
		N≥2 || throw(ErrorException("N=$N not compatible with $sym"))
		return p[1]*sin(p[2])
	elseif Base.sym_in(sym,(:z,:Z))
		N≥3 || throw(ErrorException("N=$N not compatible with $sym"))
		return p[3]
	elseif Base.sym_in(sym,(:s,))
		return p[1]
	elseif Base.sym_in(sym,(:r,))
		if N==2
			return p[1]
		elseif N==3
			return hypot(p[1],p[3])
		end
	elseif Base.sym_in(sym,(:ϕ,:φ,:𝝋,:phi))
		2≤N≤3 || throw(ErrorException("$sym only defined for 2≤N≤3, but N=$N"))
		return p[2]
	elseif Base.sym_in(sym,(:θ,:theta,:ϑ))
		N==3 || throw(ErrorException("$sym only defined for N=3, but N=$N"))
		return acos(p[3]/hypot(p[1],p[3]))
	else
		return getfield(p,sym)
	end
end
function Base.getproperty(p::Point{N,Spherical},sym::Symbol) where N
	if sym == :vec
		return getfield(p,:vector)
	elseif Base.sym_in(sym,(:x,:X))
		N≥1 || throw(ErrorException("N=$N not compatible with $sym"))
		return p[1]*sin(p[3])*cos(p[2])
	elseif Base.sym_in(sym,(:y,:Y))
		N≥2 || throw(ErrorException("N=$N not compatible with $sym"))
		return p[1]*sin(p[3])*sin(p[2])
	elseif Base.sym_in(sym,(:z,:Z))
		N≥3 || throw(ErrorException("N=$N not compatible with $sym"))
		return p[1]*cos(p[3])
	elseif Base.sym_in(sym,(:s,))
		return p[1]*sin(p[3])
	elseif Base.sym_in(sym,(:r,))
		return p[1]
	elseif Base.sym_in(sym,(:ϕ,:φ,:𝝋,:phi))
		2≤N≤3 || throw(ErrorException("$sym only defined for 2≤N≤3, but N=$N"))
		return p[2]
	elseif Base.sym_in(sym,(:θ,:theta,:ϑ))
		N==3 || throw(ErrorException("$sym only defined for N=3, but N=$N"))
		return p[3]
	else
		return getfield(p,sym)
	end
end
#1D
function Base.propertynames(::Point{1},private=false)
	if private
		return fieldnames(Point)
	else
		return (:x, :vector)
	end
end
#2D
function Base.propertynames(::Point{2,C},private=false) where C
	if private
		return fieldnames(Point)
	elseif C<:Cartesian
		return (:x, :y, :vector, :r, :ϕ)
	elseif C<:Polar
		return (:s, :ϕ, :vector, :x, :y, :r)
	end
end
#3D
function Base.propertynames(::Point{3,C},private=false) where C
	if private
		return fieldnames(Point)
	elseif C<:Cartesian
		return (:x, :y, :z, :vector, :r, :ϕ, :θ)
	elseif C<:Polar
		return (:s, :ϕ, :z, :vector, :x, :y, :r, :θ)
	elseif C<:Spherical
		return (:r, :ϕ, :θ, :vector, :x, :y, :z)
	end
end

################################################################################
# Pretty Printing

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK

# 1d (Cartesian)
function Base.show(io::IO,p::Point{1})
	print(io,"1D ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": ")
	printstyled(io, fmt("3.3f",p.x),color=PRINTED_COLOR_NUMBER)
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
	print(io, ", ϕ=")
	printstyled(io, fmt("3.3f",180*p.ϕ/π),"°",color=PRINTED_COLOR_NUMBER)
	print(io, ", θ=")
	printstyled(io, fmt("3.3f",180*p.θ/π),"°",color=PRINTED_COLOR_NUMBER)
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

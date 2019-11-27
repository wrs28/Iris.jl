"""
	module Points
"""
module Points

export Point

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK
using LinearAlgebra
using StaticArrays

struct Point{N}
	vector::SVector{N,Float64}
	Point(vec::SVector{N,T}) where {N,T<:Real} = new{N}(vec)
end

Point(x::Real,y...) = Point((x,y...))
Point(vec::NTuple{N}) where N = Point(SVector{N,Float64}(vec))
Point(p::Point) = p

"""
	struct Point{N}

`N`-Dimensional point.
Can be accessed with fields `:x`, `:y`, `:z`, `:r`, `:ϕ`, `:θ` and their variants.

------------
	Point(x,y...) -> point

"""
Point

Base.:*(p::Point,a::Number) = a*p
Base.:*(a::Number,p::Point) = Point(a*p.vec)
Base.:/(p::Point,a::Number) = Point(p.vec/a)
Base.:\(a::Number,p::Point) = p/a
Base.:+(p1::Point,p2::Point) = Point(p1.vec+p2.vec)
Base.:-(p1::Point,p2::Point) = Point(p1.vec-p2.vec)
Base.:+(p1::Point{1},p2::Number) = Point(p1.vec .+ p2)
Base.:+(p2::Number,p1::Point{1}) = Point(p1,p2)
Base.:-(p1::Point{1},p2::Number) = Point(p1.vec .- p2)
Base.:-(p2::Number,p1::Point{1}) = Point(p2 .- p1.vec)

LinearAlgebra.norm(p::Point) = norm(p.vec)

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

function Base.getproperty(p::Point{N},sym::Symbol) where N
	if sym == :vec
		return getfield(p,:vector)
	elseif Base.sym_in(sym,(:x,:X))
		N≥1 || throw(ErrorException("N=$N not compatible with x"))
		return p[1]
	elseif Base.sym_in(sym,(:y,:Y))
		N≥2 || throw(ErrorException("N=$N not compatible with y"))
		return p[2]
	elseif Base.sym_in(sym,(:z,:Z))
		N≥3 || throw(ErrorException("N=$N not compatible with z"))
		return p[3]
	elseif Base.sym_in(sym,(:r,))
		return hypot(getfield(p,:vector)...)
	elseif Base.sym_in(sym,(:ϕ,:φ,:𝝋,:phi))
		2≤N≤3 || throw(ErrorException("ϕ only defined for 2≤N≤3, but N=$N"))
		return atan(p[2],p[1])
	elseif Base.sym_in(sym,(:θ,:theta,:ϑ))
		N==3 || throw(ErrorException("θ only defined for N=3, but N=$N"))
		return acos(p[3]/hypot(getfield(p,:vector)...))
	else
		return getfield(p,sym)
	end
end

function Base.propertynames(::Point{1},private=false)
	if private
		return fieldnames(Point)
	else
		return (:x, :vector)
	end
end

################################################################################
# Pretty Printing

function Base.show(io::IO,p::Point{1})
	print(io,"1D ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": ")
	printstyled(io, p.x,color=PRINTED_COLOR_NUMBER)
end

function Base.show(io::IO,p::Point{2})
	print(io,"2D ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": x=")
	printstyled(io, p.x,color=PRINTED_COLOR_NUMBER)
	print(io, ", y=")
	printstyled(io, p.y,color=PRINTED_COLOR_NUMBER)
end

function Base.show(io::IO,p::Point{3})
	print(io,"3D ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
	print(io, ": x=")
	printstyled(io, p.x,color=PRINTED_COLOR_NUMBER)
	print(io, ", y=")
	printstyled(io, p.y,color=PRINTED_COLOR_NUMBER)
	print(io, ", z=")
	printstyled(io, p.z,color=PRINTED_COLOR_NUMBER)
end

function Base.show(io::IO,p::Point{N}) where N
	print(io, "$(N)D ")
	printstyled(io,"Point",color=PRINTED_COLOR_DARK)
end

end # module

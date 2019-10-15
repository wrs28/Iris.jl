"""
	module Points
"""
module Points

export Point

using StaticArrays

struct Point{N}
	vec::SVector{N,Float64}
	Point(vec::SVector{N,T}) where {N,T<:Real} = new{N}(vec)
end

Point(x::Real,y...) = Point((x,y...))
Point(vec::NTuple{N}) where N = Point(SVector{N,Float64}(vec))
Point(p::Point) = p
"""
	struct Point{N}

---------
	Point(x,y...) -> point

"""
Point

Base.:*(p::Point,a::Number) = a*p
Base.:*(a::Number,p::Point) = Point(a*p.vec)
Base.:+(p1::Point,p2::Point) = Point(p1.vec+p2.vec)
Base.:-(p1::Point,p2::Point) = Point(p1.vec-p2.vec)

function Base.getproperty(p::Point{N},sym::Symbol) where N
	if Base.sym_in(sym,(:x,:X))
		N≥1 || throw(ErrorException("N=$N not compatible with x"))
		return p[1]
	elseif Base.sym_in(sym,(:y,:Y))
		N≥2 || throw(ErrorException("N=$N not compatible with y"))
		return p[2]
	elseif Base.sym_in(sym,(:z,:Z))
		N≥3 || throw(ErrorException("N=$N not compatible with z"))
		return p[3]
	elseif Base.sym_in(sym,(:r,))
		return hypot(getfield(p,:vec)...)
	elseif Base.sym_in(sym,(:ϕ,:φ,:𝝋,:phi))
		2≤N≤3 || throw(ErrorException("ϕ only defined for 2≤N≤3, but N=$N"))
		return atan(p[2],p[1])
	elseif Base.sym_in(sym,(:θ,:theta,:ϑ))
		N==3 || throw(ErrorException("θ only defined for N=3, but N=$N"))
		return acos(p[3]/hypot(getfield(p,:vec)...))
	else
		return getfield(p,sym)
	end
end

Base.getindex(p::Point,index::Integer) = getfield(p,:vec)[index]

Base.ndims(::Point{N}) where N = N

Base.length(p::Point) = length(p.vec)
Base.eltype(p::Point) = Float64
Base.iterate(p::Point) = (p.vec[1],1)
Base.iterate(p::Point,state) = state ≥ length(p) ? nothing : (p.vec[state+1],state+1)

function Base.show(io::IO,p::Point{N}) where N
	print(io,N,"-dimensional point:")
	if N≤3
		N≥1 ? print(io," (x=",p.x) : nothing
		N≥2 ? print(io,", y=",p.y) : nothing
		N≥3 ? print(io,", z=",p.z) : nothing
		print(io,")")
	else
		print(io," ",p.vec)
	end
end

end # module

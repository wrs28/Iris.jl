"""
    module Lattice

for defining discrete Cartesian and Polar lattice grids used in Iris.
"""
module Lattices

using ...Defaults
using Formatting
using RecipesBase
using StaticArrays

export Lattice,
Cartesian,
Polar,
get_lattice_index,
LatticeError


struct Lattice
    type::Symbol

    dx::Float64
    dy::Float64
    θ::Float64
    x0::Float64
    y0::Float64

    dr::Float64
    dθ::Float64
    r0::Float64

    sinθ::Float64 # to avoid inefficient trig calls later
    cosθ::Float64 # to avoid inefficient trig calls later
    v1::SArray{Tuple{2},Float64,1,2}
    v2::SArray{Tuple{2},Float64,1,2}
    origin::SArray{Tuple{2},Float64,1,2}

    function Lattice(type,args...)
        if type==:Cartesian
            CartesianLat(args...)
        else
            PolarLat(args...)
        end
    end

    function CartesianLat(dx::Number,dy::Number,θ::Number=0,x0::Number=0,y0::Number=0)
        sinθ, cosθ = sincos(θ)
        R = @SMatrix [ cos(θ) -sin(θ); sin(θ)  cos(θ)]
        v1 = R*SVector{2}([1,0])
        v2 = R*SVector{2}([0,1])
        origin = SVector{2}([x0,y0])
        return new(:Cartesian,dx, dy, θ, x0, y0, NaN, NaN, NaN, sinθ, cosθ, v1, v2, origin)
    end

    function PolarLat(dr::Number,dθ::Number,θ::Number=0,x0::Number=0,y0::Number=0,r0::Number=dr/2)
        return new(:Polar,NaN,NaN,θ,x0,y0,dr,dθ,r0)
    end

    # pretty printing of Lattice
    function Base.show(io::IO, lat::Lattice)
        if lat.type==:Cartesian
            if !get(io, :sub, false)
                print(io, "Cartesian Lattice", ": \n",
                "\tdx=", fmt("3.2f",lat.dx), ", dy=", fmt("3.2f",lat.dy), ", ∠", fmt("3.2f",(mod2pi(lat.θ+π/2))*180/π), "°\n",
                "\torigin: (", fmt("3.2f",lat.x0), ", ", fmt("3.2f",lat.y0), ")")
            elseif !get(io, :sub1, false)
                print(io,
                "\n\tdx=", fmt("3.2f",lat.dx), ", dy=", fmt("3.2f",lat.dy), ", ∠", fmt("3.2f",(mod2pi(lat.θ+π/2))*180/π), "°\n",
                "\torigin: (", fmt("3.2f",lat.x0), ", ", fmt("3.2f",lat.y0), ")")
            else
                print(io,
                "\n\t\tdx=", fmt("3.2f",lat.dx),", dy=", fmt("3.2f",lat.dy), ", ∠", fmt("3.2f",(mod2pi(lat.θ+π))*180/π), "°\n",
                "\t\torigin: (", fmt("3.2f",lat.x0), ", ", fmt("3.2f",lat.y0), ")")
            end
        elseif lat.type==:Polar
            if !get(io, :sub, false)
                print(io, "Polar Lattice", ": \n",
                "\tdr=", fmt("3.2f",lat.dr), ", dθ=", fmt("3.2f",lat.dθ), "=π/", fmt("3.1f",π/lat.dθ),", ∠", fmt("3.2f",(mod2pi(lat.θ))*180/π), "°\n",
                "\torigin: (", fmt("3.2f",lat.x0), ", ", fmt("3.2f",lat.y0), ")")
            elseif !get(io, :sub1, false)
                print(io,
                "\n\tdr=", fmt("3.2f",lat.dr), ", dθ=", fmt("3.2f",lat.dθ), "=π/", fmt("3.1f",π/lat.dθ),", ∠", fmt("3.2f",(mod2pi(lat.θ))*180/π), "°\n",
                "\torigin: (", fmt("3.2f",lat.x0), ", ", fmt("3.2f",lat.y0), ")")
            else
                print(io,
                "\n\t\tdr=", fmt("3.2f",lat.dr), ", dθ=", fmt("3.2f",lat.dθ), "=π/", fmt("3.1f",π/lat.dθ),", ∠", fmt("3.2f",(mod2pi(lat.θ))*180/π), "°\n",
                "\t\torigin: (", fmt("3.2f",lat.x0), ", ", fmt("3.2f",lat.y0), ")")
            end
        else
            println(io, "unknown lattice type")
        end
    end
end
"""
    Lattices.Cartesian(dx,dy,θ=0,x0=0,y0=0) -> lattice

`x0,y0` is the origin, `θ` is the rotation

----
    (::Cartesian)(i,j) -> x,y

`x,y` are coordinates of the `i,j` lattice site (`i,j` continuous)
--------
    (lat::Cartesian)(dr=lat.dr,...)
"""
Cartesian(args...) = Lattice(:Cartesian,args...)


"""
    Lattices.Polar(dr,dθ,θ=0,r0=dr/2,x0=0,y0=0) -> lattice

`r0` is the minimum radius, `x0,y0` is the origin, `θ` is the zero angle

----
    (::Polar)(i,j) -> x,y

`x,y` are coordinates of the `i,j` lattice site (`i,j` continuous)

--------
    (lat::Polar)(dr=lat.dr,...)
"""
Polar(args...) = Lattice(:Polar,args...)


"""
    (::Lattice)(i,j) -> x,y
"""
function (lat::Lattice)(i::Real,j::Real)
    if lat.type==:Cartesian
        v = lat.origin + i*lat.v1*lat.dx + j*lat.v2*lat.dy
        return v[1],v[2]
    elseif lat.type==:Polar
        r = lat.r0 + i*lat.dr
        θ = lat.θ  + j*lat.dθ
        v = r.*sincos(θ)
        return lat.x0+v[2],lat.y0+v[1]
    end
end

# functor wrappers for both Cartesian and Polar lattices
(lat::Lattice)(ijk::Tuple) = lat(ijk...)
(lat::Lattice)(k::CartesianIndex) = lat(k[1],k[2])
function (lat::Lattice)(i::Real,j::Real,v::AbstractArray)
    v[:] = lat(i,j)
    return nothing
end
function (lat::Lattice)(x::AbstractArray,y::AbstractArray)
    X = Array{Float64}(undef,size(x))
    Y = Array{Float64}(undef,size(y))
    for i ∈ eachindex(X)
        X[i], Y[i] = lat(x[i],y[i])
    end
    return X,Y
end

function (lat::Lattice)(;dx=lat.dx, dy=lat.dy, θ=lat.θ, x0=lat.x0, y0=lat.y0, dr=lat.dr, dθ=lat.dθ, r0=dr/2)
    if lat.type==:Cartesian
        return Cartesian(dx,dy,θ,x0,y0)
    else
        return Polar(dr,dθ,θ,x0,y0,r0)
    end
end


struct LatticeError <: Exception
    lat::Lattice
end
function Base.showerror(io::IO,err::LatticeError)
    print(io,"LatticeError: unrecognized lattice type ")
    print(io, err.lat.type)
    print(io, ", must be one of :Cartesian or :Polar")
end

"""
    get_lattice_index(lattice,x,y) -> nx,ny

returns indices `nx`,`ny` such that `(x,y) = (x0,y0) + nx*v1*dx + ny*v2*dx`
"""
function get_lattice_index(lat::Lattice,x::AbstractArray,y::AbstractArray)
    @assert size(x)==size(y) "x and y must have the same size"
    X = Array{Float64}(undef,size(x)...)
    Y = Array{Float64}(undef,size(y)...)
    for k ∈ eachindex(X)
        X[k],Y[k] = get_lattice_index(lat,x[k],y[k])
    end
    return X,Y
end
function get_lattice_index(lat::Lattice,x::Real,y::Real)
    i,j = lat.type==:Cartesian ? get_lattice_index_cartesian(lat,x,y) : get_lattice_index_polar(lat,x,y)
    return i,j
end
function get_lattice_index_cartesian(lat::Lattice,x::Real,y::Real)
    dx, dy, v1, v2, x0, y0 = lat.dx, lat.dy, lat.v1, lat.v2, lat.x0, lat.y0
    x += -x0
    y += -y0
    rv1 = v1[1]*x + v1[2]*y
    rv2 = v2[1]*x + v2[2]*y
    return rv1/dx, rv2/dy
end
function get_lattice_index_polar(lat::Lattice,x::Real,y::Real)
    dr, dθ, x0, y0 = lat.dr, lat.dθ, lat.x0, lat.y0
    x += -x0
    y += -y0
    r, θ = hypot(y,x), atan(y,x)+π
    return (r-lat.r0)/dr, (θ-lat.θ)/dθ
end


"""
    plot(lattice, N=[4,4])

scatter plot of lattice vectors from -`N[i]`:`+N[i]` for each dimension i.

for polar lattices, plots all angles out to radial site `N[1]`.

If `N` is a scalar, uses same span for both dimensions.
"""
@recipe f(lattice::Lattice, N::Int=4) = (lattice, [N,N])
@recipe function f(lattice::Lattice, N::Array{Int,1})
    lattice.type==:Polar ? N[2]=floor(Int,2π/lattice.dθ-eps(2π/lattice.dθ)) : nothing
    seriestype --> :scatter
    legend --> false
    aspect_ratio --> 1
    @series begin
        x, y = Float64[], Float64[]
        for n ∈ 0:N[1], m ∈ 0:N[2]
            temp = lattice(n,m)
            push!(x,temp[1]); push!(y,temp[2])
        end
        (x, y)
    end
end


end # module

"""
    module Shapes

Definitions of various shapes.

All shapes have the signature: `shape(parameters::Tuple, x0, y0, θ=0)`, and some
have keyword argument `reference` which determines which point `(x0,y0)` is referring to.

The length of `parameters` depends on the shape. For example, `Circle` has just one
parameter, the radius, while `Ellipse` has two, the semi-major and semi-minor axes.
"""
module Shapes

export Circle,
Ellipse,
DeformedDisk,
Annulus,
Square,
Rectangle,
Parallelogram,
normal_distance,
get_number_of_sides,
rotate,
unrotate

export AbstractShape,
AbstractDisk,
AbstractQuadrilateral

using ...Defaults
using Formatting
using LinearAlgebra
using NLopt
using RecipesBase
using StaticArrays

abstract type AbstractShape{N} end # N is the number of sides
abstract type AbstractDisk <: AbstractShape{1} end
abstract type AbstractQuadrilateral <: AbstractShape{4} end

get_number_of_sides(s::AbstractShape{N}) where N = N

"""
    struct Circle <: AbstractDisk <: AbstractShape{1}

    Circle((R,), x0, y0, θ=0; reference=:center) -> circle
    Circle(R, x0, y0, θ=0; reference=:center) -> circle

`R` radius
`x0, y0` is location of `reference`

-------------------

    (::Circle)(x,y) --> is_in_circle::Bool
"""
struct Circle <: AbstractDisk
    R::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    sinθ::Float64
    cosθ::Float64
    corners::Tuple{}
    models::Array{Opt,1}

    Circle((R,), args...; kwargs...) = Circle(R[1], args...; kwargs...)
    function Circle(R::Real,x0::Real,y0::Real,θ::Real=0; reference::Symbol=:center)
        if reference ∈ [:bottom, :Bottom, :b, :B]
            x0,y0 = (x0,y0) .+ rotate( 0, R,cos(θ),sin(θ))
        elseif reference ∈ [:top, :Top, :T, :t]
            x0,y0 = (x0,y0) .+ rotate( 0,-R,cos(θ),sin(θ))
        elseif reference ∈ [:left, :Left, :L, :l]
            x0,y0 = (x0,y0) .+ rotate( R, 0,cos(θ),sin(θ))
        elseif reference ∈ [:right, :Right, :R, :r]
            x0,y0 = (x0,y0) .+ rotate(-R, 0,cos(θ),sin(θ))
        end
        c = new(R,x0,y0,θ,sin(θ),cos(θ))
        models = build_models(c)
        return new(R,x0,y0,θ,sin(θ),cos(θ),(),models)
    end

    (c::Circle)(x,y) = hypot((x-c.x0), (y-c.y0)) < c.R-eps(c.R)

    function Base.show(io::IO, circle::Circle)
        print(io, "Circle(R=", fmt("2.2f",circle.R), ", x0=", fmt("2.2f",circle.x0), ", y0=", fmt("2.2f",circle.y0), ")")
    end
end
"""
    perimeter(t,shape,side) -> x,y,dx/dt,dy/dt
"""
function perimeter(t::Number,c::Circle,side::Int)
    sn,cs = sincos(t)
    x = c.x0 + c.R*cs
    y = c.y0 + c.R*sn
    ẋ = -c.R*sn
    ẏ =  c.R*cs
    return x,y,ẋ,ẏ
end


"""
    struct Ellipse <: AbstractDisk <: AbstractShape{1}

    Ellipse((a,b), x0, y0, θ=0; reference=:center) -> ellipse
    Ellipse(a, b, x0, y0, θ=0; reference=:center) -> ellipse


`a, b` axis lengths
`x0, y0` is location of `reference`
`θ` angle of `a` axis

----------------

    (::Ellipse)(x,y) -> is_in_ellipse::Bool
"""
struct Ellipse <: AbstractDisk
    a::Float64
    b::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    sinθ::Float64
    cosθ::Float64
    corners::Tuple{}
    models::Array{Opt,1}

    Ellipse((a,b),args...;kwargs...) = Ellipse(a,b,args...;kwargs...)
    function Ellipse(a::Number,b::Number,x0::Number,y0::Number,θ::Number=0; reference::Symbol=:center)
        if reference ∈ [:bottom, :Bottom, :b, :B]
            x0,y0 = (x0,y0) .+ rotate( 0, b,cos(θ),sin(θ))
        elseif reference ∈ [:top, :Top, :T, :t]
            x0,y0 = (x0,y0) .+ rotate( 0,-b,cos(θ),sin(θ))
        elseif reference ∈ [:left, :Left, :L, :l]
            x0,y0 = (x0,y0) .+ rotate( a, 0,cos(θ),sin(θ))
        elseif reference ∈ [:right, :Right, :R, :r]
            x0,y0 = (x0,y0) .+ rotate(-a, 0,cos(θ),sin(θ))
        end
         e = new(a,b,x0,y0,θ,sin(θ),cos(θ))
         models = build_models(e)
         return new(a,b,x0,y0,θ,sin(θ),cos(θ),(),models)
    end

    Ellipse(c::Circle) = Ellipse(c.R,c.R,c.x0,c.y0,c.θ)

    function (e::Ellipse)(x,y)
        xurot, yurot = unrotate(x, y, e)
        return hypot((xurot-e.x0)/e.a, (yurot-e.y0)/e.b) < 1-eps()
    end

    function Base.show(io::IO, ellipse::Ellipse)
        print(io, "Ellipse(a=", fmt("2.2f",ellipse.a), ", b=", fmt("2.2f",ellipse.b), ", x0=", fmt("2.2f",ellipse.x0), ", y0=", fmt("2.2f",ellipse.y0), ", θ=∠", fmt("3.2f",(mod2pi(ellipse.θ))*180/π), "°)")
    end
end
function perimeter(t::Number,e::Ellipse,side::Int)
    sn,cs = sincos(t)
    xur,yur =    e.a*cs, e.b*sn
    xpur,ypur = -e.a*sn, e.b*cs
    x,y = (e.x0,e.y0) .+ rotate(xur,yur,e.cosθ,e.sinθ)
    ẋ,ẏ = rotate(xpur,ypur,e.cosθ,e.sinθ)
    return x,y,ẋ,ẏ
end


"""
    struct DeformedDisk{N} <: AbstractDisk <: AbstractShape{1}

    DeformedDisk{N}(R, x0, y0, M, a, φ, θ=0) -> deformeddisk
    DeformedDisk{N}((R,M,a,φ), x0, y0, θ=0) -> deformeddisk

`R` is radius,
`x0` and `y0` is center of circle of radius `R`,
`M` is array of length `N` of multipole integers,
`a` is array of length `N` amplitudes,
`φ` is array of length `N` of angles
`θ` is overall rotation angle

-----------------

    (::DeformedDisk)(x,y) -> is_in_disk::Bool
"""
struct DeformedDisk{N} <: AbstractDisk
    R::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    M::SArray{Tuple{N},Int,1,N}
    a::SArray{Tuple{N},Float64,1,N}
    φ::SArray{Tuple{N},Float64,1,N}
    sinθ::Float64
    cosθ::Float64
    corners::Tuple{}
    models::Array{Opt,1}

    DeformedDisk{n}((R,M,a,φ),x0,y0,args...;kwargs...) where n= DeformedDisk{n}(R,x0,y0,M,a,φ,args...;kwargs...)
    function DeformedDisk{n}(R::Number,x0::Number,y0::Number,M,a,φ,θ=0) where n
        @assert n==length(M)==length(a)==length(φ) "parameter N in DeformedDisk{N}(...) must be equal to length(M)"
        d = new{n}(R,x0,y0,θ,SVector{n}(M),SVector{n}(a),SVector{n}(φ.+θ),sin(θ),cos(θ))
        models = build_models(d)
        return new{n}(R,x0,y0,d.θ,d.M,d.a,d.φ,d.sinθ,d.cosθ,(),models)
    end

    DeformedDisk(c::Circle) = DeformedDisk{0}(c.R,c.x0,c.y0,[],[],[],c.θ)

    function (d::DeformedDisk)(x,y)
        θ = atan(y-d.y0,x-d.x0)
        r = hypot(x-d.x0,y-d.y0) < d.R + sum(d.a.*(map((m,φ)->cos(m*(θ-φ)),d.M,d.φ))) - eps(d.R)
        return r
    end

    function Base.show(io::IO, d::DeformedDisk{N}) where N
         print(io,"DeformedDisk with $N multipoles")
    end
end
function perimeter(t::Number,d::DeformedDisk,side::Int)
    r = d.R + sum(d.a.*(map((m,φ)->cos(m*(t-φ)),d.M,d.φ)))
    ṙ = sum(d.a.*(map((m,φ)->-m*sin(m*(t-φ)),d.M,d.φ)))
    sn,cs = sincos(t)
    x,y = d.x0 + r*cs, d.y0 + r*sn
    ẋ = ṙ*cs - r*sn
    ẏ = ṙ*sn + r*cs
    return x,y,ẋ,ẏ
end

# plotting for all disks
@recipe function f(d::AbstractDisk)
    θ = LinRange(0,2π,201)
    x = Array{Float64,1}(undef,length(θ))
    y = Array{Float64,1}(undef,length(θ))
    for i ∈ eachindex(θ)
        x[i],y[i] = perimeter(θ[i],d,1)
    end
    alpha --> 0
    seriestype --> :path
    fillcolor --> SHAPE_COLOR
    fillrange --> d.y0
    fillalpha --> SHAPE_FILL_ALPHA
    aspect_ratio --> 1
    legend-->false
    (x,y)
end


"""
    struct Annulus <: AbstractShape{2}

    Annulus((R1,R2),x0,y0,θ=0) -> annulus
    Annulus(R1, R2, x0, y0, θ=0) -> annulus

`R1` is inner radius, `R2` outer
`x0, y0` is center

-------------
    (::Annulus)(x,y) -> is_in_annulus::Bool
"""
struct Annulus <: AbstractShape{2}
    R1::Float64
    R2::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    sinθ::Float64
    cosθ::Float64
    corners::Tuple{}
    models::Array{Opt,1}

    Annulus((R1,R2),args...;kwargs...) = Annulus(R1,R2,args...;kwargs...)
    function Annulus(R1::Number,R2::Number,x0::Number,y0::Number,θ::Number=0)
        @assert R1≤R2 "R1=$R1, R2=$R2 must satisfy R1≤R2"
        a = new(R1,R2,x0,y0,0,0,1)
        models = build_models(a)
        return new(R1,R2,x0,y0,0,sin(θ),cos(θ),(),models)
    end

    (c::Annulus)(x,y) = c.R1+eps(c.R1) < hypot( (x-c.x0),(y-c.y0) ) < c.R2-eps(c.R2)

    function Base.show(io::IO, a::Annulus)
        print(io, "Annulus(R1=$(fmt("2.2f",a.R1)), R2=$(fmt("2.2f",a.R2)), x0=$(fmt("2.2f",a.x0)), y0=$(fmt("2.2f",a.y0)))")
    end
end
function perimeter(t::Number,a::Annulus,side::Int)
    sn,cs = sincos(t)
    R = side==1 ? a.R1 : a.R2
    x = a.x0 + R*cs
    y = a.y0 + R*sn
    ẋ = -R*sn
    ẏ =  R*cs
    return x,y,ẋ,ẏ
end
@recipe function f(cr::Annulus)
    @series Circle(cr.R1,cr.x0,cr.y0)
    @series Circle(cr.R2,cr.x0,cr.y0)
end


"""
    struct Square <: AbstractQuadrilateral <: AbstractShape{4}

    Square((a,),x0,y0,θ=0; reference=:center) -> square
    Square(a,x0,y0,θ=0); reference=:center -> square

`a` length of sides
`x0, y0` is location of `reference`
`θ` angle of `a` side

------------

    (::Square)(x,y) -> is_in_square::Bool
"""
struct Square <: AbstractQuadrilateral
    a::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    sinθ::Float64
    cosθ::Float64
    corners::NTuple{4,NTuple{2,Float64}}
    models::Array{Opt,1}

    Square((a,),args...;kwargs...) = Square(a,args...;kwargs...)
    function Square(a::Number,x0::Number,y0::Number,θ::Number=0; reference::Symbol=:center)
        if reference ∈ [:bottom, :Bottom, :b, :B]
            x0,y0 = (x0,y0) .+ rotate( 0, a/2,cos(θ),sin(θ))
        elseif reference ∈ [:top, :Top, :T, :t]
            x0,y0 = (x0,y0) .+ rotate( 0,-a/2,cos(θ),sin(θ))
        elseif reference ∈ [:left, :Left, :L, :l]
            x0,y0 = (x0,y0) .+ rotate( a/2, 0,cos(θ),sin(θ))
        elseif reference ∈ [:right, :Right, :R, :r]
            x0,y0 = (x0,y0) .+ rotate(-a/2, 0,cos(θ),sin(θ))
        elseif reference ∈ [:topright, :TopRight, :TR, :tr, :ne, :NE]
            x0,y0 = (x0,y0) .+ rotate(-a/2,-a/2,cos(θ),sin(θ))
        elseif reference ∈ [:topleft, :TopLeft, :TL, :tl, :nw, :NW]
            x0,y0 = (x0,y0) .+ rotate( a/2,-a/2,cos(θ),sin(θ))
        elseif reference ∈ [:bottomleft, :BottomLeft, :BL, :bl, :sw, :SW]
            x0,y0 = (x0,y0) .+ rotate( a/2, a/2,cos(θ),sin(θ))
        elseif reference ∈ [:bottomright, :BottomRight, :BR, :br, :se, :SE]
            x0,y0 = (x0,y0) .+ rotate(-a/2, a/2,cos(θ),sin(θ))
        end
        s = new(a,x0,y0,θ,sin(θ),cos(θ))
        c1 = rotate(s.x0+s.a/2,s.y0+s.a/2, s)
        c2 = rotate(s.x0-s.a/2,s.y0+s.a/2, s)
        c3 = rotate(s.x0-s.a/2,s.y0-s.a/2, s)
        c4 = rotate(s.x0+s.a/2,s.y0-s.a/2, s)
        models = build_models(s)
        return new(a,x0,y0,θ,sin(θ),cos(θ),(c1,c2,c3,c4),models)
    end

    function (s::Square)(x,y)
        xurot, yurot = unrotate(x, y, s)
        return -s.a/2+10eps(s.a/2) < xurot-s.x0 < s.a/2-10eps(s.a/2)  &&  -s.a/2+10eps(s.a/2) < yurot-s.y0 < s.a/2-10eps(s.a/2)
    end

    function Base.show(io::IO, square::Square)
        print(io, "Square(a=$(fmt("2.2f",square.a)), x0=$(fmt("2.2f",square.x0)), y0=$(fmt("2.2f",square.y0)), θ=∠", fmt("3.2f",(mod2pi(square.θ))*180/π), "°)")
    end
end
function perimeter(t::Number,s::Square,side::Int)
    if side==1
        x = -s.a/2
        y = -s.a/2 + t*s.a
        ẋ = 0
        ẏ = s.a
    elseif side==2
        x =  s.a/2
        y = -s.a/2 + t*s.a
        ẋ = 0
        ẏ = s.a
    elseif side==3
        x = -s.a/2 + t*s.a
        y = -s.a/2
        ẋ = s.a
        ẏ = 0
    elseif side==4
        x = -s.a/2 + t*s.a
        y =  s.a/2
        ẋ = s.a
        ẏ = 0
    end
    x += s.x0
    y += s.y0
    xrot,yrot = rotate(x,y,s)
    ẋrot,ẏrot = rotate(ẋ,ẏ,s.cosθ,s.sinθ)
    return xrot,yrot,ẋrot,ẏrot
end


"""
    struct Rectangle <: AbstractQuadrilateral <: AbstractShape{4}

    Rectangle((a,b), x0, y0, θ=0; reference=:center) -> rectangle
    Rectangle(a, b, x0, y0, θ=0; reference=:center) -> rectangle

`a, b` length of sides
`x0, y0` is location of `reference`
`θ` angle of `a` side

-----------------------

    (::Rectangle)(x,y) -> is_in_rectangle::Bool
"""
struct Rectangle <: AbstractQuadrilateral
    a::Float64
    b::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    sinθ::Float64
    cosθ::Float64
    corners::NTuple{4,NTuple{2,Float64}}
    models::Array{Opt,1}

    Rectangle((a,b),args...;kwargs...) = Rectangle(a,b,args...;kwargs...)
    function Rectangle(a::Number,b::Number,x0::Number,y0::Number,θ::Number=0; reference::Symbol=:center)
        if reference ∈ [:bottom, :Bottom, :b, :B]
            x0,y0 = (x0,y0) .+ rotate( 0, b/2,cos(θ),sin(θ))
        elseif reference ∈ [:top, :Top, :T, :t]
            x0,y0 = (x0,y0) .+ rotate( 0,-b/2,cos(θ),sin(θ))
        elseif reference ∈ [:left, :Left, :L, :l]
            x0,y0 = (x0,y0) .+ rotate( a/2, 0,cos(θ),sin(θ))
        elseif reference ∈ [:right, :Right, :R, :r]
            x0,y0 = (x0,y0) .+ rotate(-a/2, 0,cos(θ),sin(θ))
        elseif reference ∈ [:topright, :TopRight, :TR, :tr, :ne, :NE]
            x0,y0 = (x0,y0) .+ rotate(-a/2,-b/2,cos(θ),sin(θ))
        elseif reference ∈ [:topleft, :TopLeft, :TL, :tl, :nw, :NW]
            x0,y0 = (x0,y0) .+ rotate( a/2,-b/2,cos(θ),sin(θ))
        elseif reference ∈ [:bottomleft, :BottomLeft, :BL, :bl, :sw, :SW]
            x0,y0 = (x0,y0) .+ rotate( a/2, b/2,cos(θ),sin(θ))
        elseif reference ∈ [:bottomright, :BottomRight, :BR, :br, :se, :SE]
            x0,y0 = (x0,y0) .+ rotate(-a/2, b/2,cos(θ),sin(θ))
        end
        r = new(a,b,x0,y0,θ,sin(θ),cos(θ))
        models = build_models(r)
        corners_x, corners_y = Array{Float64}(undef,4), Array{Float64}(undef,4)
        c1 = (r.x0,r.y0) .+ rotate( r.a/2, r.b/2,r.cosθ,r.sinθ)
        c2 = (r.x0,r.y0) .+ rotate(-r.a/2, r.b/2,r.cosθ,r.sinθ)
        c3 = (r.x0,r.y0) .+ rotate(-r.a/2,-r.b/2,r.cosθ,r.sinθ)
        c4 = (r.x0,r.y0) .+ rotate( r.a/2,-r.b/2,r.cosθ,r.sinθ)
        return new(a,b,x0,y0,θ,sin(θ),cos(θ),(c1,c2,c3,c4),models)
    end

    Rectangle(s::Square) = Rectangle(s.a,s.a,s.x0,s.y0,s.θ)

    function (r::Rectangle)(x,y)
        xurot, yurot = unrotate(x, y, r)
        return -r.a/2+20eps(r.a/2) < xurot-r.x0 < r.a/2-20eps(r.a/2)  &&  -r.b/2+20eps(r.b/2) < yurot-r.y0 < r.b/2-20eps(r.b/2)
    end

    function Base.show(io::IO, rect::Rectangle)
        print(io, "Rectangle(a=$(fmt("2.2f",rect.a)), b=$(fmt("2.2f",rect.b)), x0=$(fmt("2.2f",rect.x0)), y0=$(fmt("2.2f",rect.y0)), θ=∠", fmt("3.2f",(mod2pi(rect.θ))*180/π), "°)")
    end
end
function perimeter(t::Number,s::Rectangle,side::Int)
    if side==1
        x = -s.a/2
        y = -s.b/2 + t*s.b
        ẋ = 0
        ẏ = s.b
    elseif side==2
        x =  s.a/2
        y = -s.b/2 + t*s.b
        ẋ = 0
        ẏ = s.b
    elseif side==3
        x = -s.a/2 + t*s.a
        y = -s.b/2
        ẋ = s.a
        ẏ = 0
    elseif side==4
        x = -s.a/2 + t*s.a
        y =  s.b/2
        ẋ = s.a
        ẏ = 0
    end
    x += s.x0
    y += s.y0
    xrot,yrot = rotate(x,y,s)
    ẋrot,ẏrot = rotate(ẋ,ẏ,s.cosθ,s.sinθ)
    return xrot,yrot,ẋrot,ẏrot
end


# """
#     Parallelogram(a, b, α, x0, y0, θ)
# """
# struct Parallelogram <: AbstractParallelogram
#     a::Float64
#     b::Float64
#     α::Float64
#     x0::Float64
#     y0::Float64
#     θ::Float64
#     sinθ::Float64
#     cosθ::Float64
#     tanα::Float64
#     cosα::Float64
#     models::Array{Opt,1}
#
#     function Parallelogram(a::Number, b::Number, α::Number, x0::Number, y0::Number, θ::Number)
#         new(a,b,α,x0,y0,θ,cos(θ),sin(θ),tan(α),cos(α))
#     end
#
#     function (p::Parallelogram)(x,y)
#         xrot, yrot = rotate(x-p.x0, y-p.y0, p.cosθ, p.sinθ)
#         return p.tanα*(xrot-p.a) < yrot < p.tanα*xrot  &&  0 < yrot < p.cosα*p.b
#     end
#
#     Base.show(io::IO, par::Parallelogram) = print(io, "Parallelogram(a=$(fmt("2.2f",par.a)), b=$(fmt("2.2f",par.b)), α=$(fmt("2.2f",par.α)), x0=$(fmt("2.2f",par.x0)), y0=$(fmt("2.2f",par.y0)), θ=$(fmt("2.2f",par.θ)))")
#
#     @recipe function f(pg::Parallelogram)
#         x = cumsum([0, pg.a, +pg.b*cos(pg.α), -pg.a, -pg.b*cos(pg.α)])
#         y = cumsum([0, 0, +pg.b*sin(pg.α), 0, -pg.b*sin(pg.α)])
#         alpha --> 0
#         seriestype --> :path
#         fillcolor --> SHAPE_COLOR
#         fillrange --> 0
#         fillalpha --> SHAPE_FILL_ALPHA
#         aspect_ratio --> 1
#         legend --> false
#         xrot, yrot = rotate(x,y,pg.cosθ,-pg.sinθ)
#         pg.x0 .+ xrot, pg.y0 .+ yrot
#     end
# end

# plot all quadrilaterals
@recipe function f(q::AbstractQuadrilateral)
    t = LinRange(0,1,201)
    x1,y1 = Array{Float64}(undef,length(t)), Array{Float64}(undef,length(t))
    x2,y2 = Array{Float64}(undef,length(t)), Array{Float64}(undef,length(t))
    x3,y3 = Array{Float64}(undef,length(t)), Array{Float64}(undef,length(t))
    x4, y4 = Array{Float64}(undef,length(t)), Array{Float64}(undef,length(t))
    for i ∈ eachindex(t)
        x1[i],y1[i] = perimeter(t[i],q,1)
        x2[i],y2[i] = perimeter(t[i],q,2)
        x3[i],y3[i] = perimeter(t[i],q,3)
        x4[i],y4[i] = perimeter(t[i],q,4)
    end
    alpha --> 0
    seriestype --> :path
    fillcolor --> SHAPE_COLOR
    fillrange --> q.y0
    fillalpha --> SHAPE_FILL_ALPHA
    aspect_ratio --> 1
    legend --> false
    vcat(x1,x4,reverse(x3)),vcat(y1,y4,reverse(y3))
end


"""
    rotate(x,y,cosθ,sinθ) -> xrot, yrot

rotate by `θ` about origin

----------------
    rotate(x,y,shape) -> xrot, yrot

rotate about `shape`'s center by `shape`'s angle
"""
rotate(x::Real,y::Real,cosθ::Real,sinθ::Real) = iszero(sinθ) ? (x,y) : (cosθ*x-sinθ*y, sinθ*x+cosθ*y)
rotate(x::Real,y::Real,s::AbstractShape) = (s.x0,s.y0) .+ rotate(x-s.x0,y-s.y0,s.cosθ,s.sinθ)
function rotate(x::AbstractArray,y::AbstractArray,cosθ::Real,sinθ::Real)
    xrot = Array{Float64}(undef,size(xyrot))
    yrot = Array{Float64}(undef,size(xyrot))
    rotate!(xrot,yrot,x,y,cosθ,sinθ)
    return xrot, yrot
end
function rotate(x::AbstractArray,y::AbstractArray,s::AbstractShape)
    xrot,yrot = rotate(x.-s.x0,y.-y.s0,s.cosθ,s.sinθ)
    for i ∈ eachindex(xrot)
        xrot[i] += s.x0
        yrot[i] += s.y0
    end
    return xrot,yrot
end
function rotate!(xrot::AbstractArray,yrot::AbstractArray,x::AbstractArray,y::AbstractArray,cosθ::Real,sinθ::Real)
    for i ∈ eachindex(x)
        xrot[i],yrot[i] = rotate(x[i],y[i],cosθ,sinθ)
    end
    return nothing
end


"""
    unrotate(x,y,cosθ,sinθ) -> xrot, yrot

rotate by `-θ` about origin

----------------
    unrotate(x,y,shape) -> xrot, yrot

rotate about `shape`'s center by the negative of `shape`'s angle
"""
unrotate(x,y,cosθ::Real,sinθ::Real) = rotate(x,y,cosθ,-sinθ)
unrotate(x::Real,y::Real,s::AbstractShape) = (s.x0,s.y0) .+ unrotate(x-s.x0,y-s.y0,s.cosθ,s.sinθ)
function unrotate(x::AbstractArray,y::AbstractArray,s::AbstractShape)
    xrot,yrot = unrotate(x.-s.x0,y.-y.s0,s.cosθ,s.sinθ)
    for i ∈ eachindex(xrot)
        xrot[i] += s.x0
        yrot[i] += s.y0
    end
    return xrot,yrot
end


function build_models(s::AbstractShape{N}) where N
    models = Array{Opt}(undef,N)
    for i ∈ eachindex(models)
        models[i] = Opt(NL_NORMAL_ALGORITHM, 1)
        models[i].xtol_rel = NL_NORMAL_XTOL
        models[i].ftol_rel = NL_NORMAL_FTOL
        models[i].maxeval = NL_NORMAL_MAXEVAL
        if typeof(s)<:AbstractQuadrilateral
            models[i].lower_bounds = 0
            models[i].upper_bounds = 1
        end
    end
    return models
end


function model_objective(v::Vector,grad::Vector,X::Number,Y::Number,s::AbstractShape,side::Int)
    x,y,ẋ,ẏ = perimeter(v[1],s,side)
    if length(grad) > 0
        grad[1] = 2(x-X)*ẋ + 2(y-Y)*ẏ
    end
    return (x-X)^2 + (y-Y)^2
end


"""
    normal_distance(shape,x,y) -> normal, tangent, distance, side

Get the `normal` vector from the surface of `shape` to a point `x,y`,
and the associated `tangent` vector and distance from surface `side`.

The returned arguments are sorted in ascending order by `distance`, so that, e.g.,
`normal[1]` is the normal vector associated with the closest side.
"""
function normal_distance(s::TS,x::Number,y::Number) where TS<:AbstractShape
    n = Array{Array{Float64,1}}(undef,length(s.models))
    d = Array{Float64}(undef,length(s.models))
    sides = 1:length(s.models)
    for i ∈ eachindex(s.models)
        model = s.models[i]
        model.min_objective = (v,g)->model_objective(v,g,x,y,s,i)
        minf,mint,ret = optimize(model,[t0_shape(x,y,s)])
        X,Y = perimeter(mint[1],s,i)
        n[i] = [x-X, y-Y]
        n[i] = n[i]/norm(n[i])
        d[i] = sqrt(minf)
    end
    perm = sortperm(d)
    n, d = n[perm], d[perm]
    NORMAL = Array{Array{Float64,1}}(undef,length(s.models))
    TANGENT = Array{Array{Float64,1}}(undef,length(s.models))
    for i ∈ eachindex(NORMAL)
        nx, ny = n[i][1],n[i][2]
        NORMAL[i] = [nx, ny]
        TANGENT[i] = [ny,-nx]
    end
    return NORMAL, TANGENT, d, sides[perm]
end

t0_shape(x::Number,y::Number,s::AbstractQuadrilateral) = .5
t0_shape(x::Number,y::Number,s::AbstractDisk) = atan(y-s.y0,x-s.x0)


end #module

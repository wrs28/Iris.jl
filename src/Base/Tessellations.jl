"""
    module Tessellations

a convenience wrapper to pacakge VoronoiDelaunay, with main export Tessellation.
"""
module Tessellations

export Tessellation,
tessellation_coords,
tessellation_inv_coords,
get_surrounding_sites

using ...Defaults
using RecipesBase
using VoronoiDelaunay


"""
    struct Tessellation

    Tessellation(x,y) -> t

create a Delaunay tessellation out of the points `x,y` (passed as vectors)
"""
struct Tessellation
    tessellation::DelaunayTessellation2D{Point2D}
    x::Array{Float64,1}
    y::Array{Float64,1}
    xy::Array{Point2D,1}
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64

    function Tessellation(x,y)
        xmin, xmax = extrema(x)
        xmin_new = (xmin+xmax)/2 - TESSELLATION_EXPANSION_FACTOR*(xmax-xmin)/2
        xmax_new = (xmin+xmax)/2 + TESSELLATION_EXPANSION_FACTOR*(xmax-xmin)/2
        ymin, ymax = extrema(y)
        ymin_new = (ymin+ymax)/2 - TESSELLATION_EXPANSION_FACTOR*(ymax-ymin)/2
        ymax_new = (ymin+ymax)/2 + TESSELLATION_EXPANSION_FACTOR*(ymax-ymin)/2
        return Tessellation(x,y,xmin_new,xmax_new,ymin_new,ymax_new)
    end
    function Tessellation(x,y,xmin,xmax,ymin,ymax)
        tessellation = DelaunayTessellation(length(x))
        xt, yt = tessellation_coords(x,y,xmin,xmax,ymin,ymax)
        p = Point.(xt,yt)
        push!(tessellation,deepcopy(p))
        return new(tessellation,x,y,p,xmin,xmax,ymin,ymax)
    end

    function Base.show(io::IO,t::Tessellation)
        Np = length(t.xy)
        Nt = length(t.tessellation._trigs)
        print(io,"Tessellation:")
        print(io," ",Np," points,")
        print(io," ",Nt," triangles")
    end
end


"""
    tessellation_coords(x,y,t::Tessellation) -> X,Y

map the point(s) `x,y` to the square ((1,2),(1,2)), needed for the VoronoiDelaunay pacakage.
"""
tessellation_coords(x,y,t::Tessellation) = tessellation_coords(x,y,t.xmin,t.xmax,t.ymin,t.ymax)
function tessellation_coords(x,y,args...)
    XT = Array{Float64}(undef,length(x))
    YT = Array{Float64}(undef,length(x))
    tessellation_coords!(XT,YT,x,y,args...)
    return XT, YT
end
function tessellation_coords!(XT,YT,x,y,args...)
    @inbounds for i ∈ eachindex(XT)
        XT[i],YT[i] = tessellation_coords(x[i],y[i],args...)
    end
    return nothing
end
function tessellation_coords(x::Real,y::Real,xmin,xmax,ymin,ymax)
    xt = min_coord + (max_coord-min_coord)*(x - xmin)/(xmax-xmin)
    yt = min_coord + (max_coord-min_coord)*(y - ymin)/(ymax-ymin)
    return xt, yt
end


"""
    tessellation_inv_coords(x,y,t::Tessellation) -> X,Y

undo the mapping to ((1,2),(1,2)) that is needed for the VoronoiDelaunay pacakage.
"""
tessellation_inv_coords(x,y,t::Tessellation) = tessellation_inv_coords(x,y,t.xmin,t.xmax,t.ymin,t.ymax)
function tessellation_inv_coords(xt,yt,xmin,xmax,ymin,ymax)
    X = Array{Float64}(undef,length(xt))
    Y = Array{Float64}(undef,length(xt))
    tessellation_inv_coords!(X,Y,xt,yt,xmin,xmax,ymin,ymax)
    return X, Y
end
function tessellation_inv_coords!(X,Y,xt,yt,xmin,xmax,ymin,ymax)
    @inbounds for i ∈ eachindex(X)
        X[i],Y[i] = tessellation_inv_coords(xt[i],yt[i],xmin,xmax,ymin,ymax)
    end
    return nothing
end
function tessellation_inv_coords(xt::Real,yt::Real,xmin,xmax,ymin,ymax)
    x = xmin + (xmax-xmin)*(xt - min_coord)/(max_coord-min_coord)
    y = ymin + (ymax-ymin)*(yt - min_coord)/(max_coord-min_coord)
    return x, y
end


"""
    get_surrounding_sites(t,x,y[,deep::Bool]) -> inds

`t` is a Tessellation, `x,y` is a positions.
returns indices of vertices in `t` that surround the point `x,y`

If `deep=true`, it goes 2-layers deep, meaning that it looks at the surrounding
triangles as well as those surrounding them.
"""
@inline function get_surrounding_sites(t::Tessellation,x::Real,y::Real,deep::Bool=false)
    t0 = locate(t,x,y)
    a, b, c = t0._a, t0._b, t0._c

    P = vcat(a,b,c)

    if isexternal(t0) || deep

        ta = movea(t,t0); Pa = [ta._a,ta._b,ta._c]
        tb = moveb(t,t0); Pb = [tb._a,tb._b,tb._c]
        tc = movec(t,t0); Pc = [tc._a,tc._b,tc._c]
        P = vcat(P,Pa,Pb,Pc)

        if deep
            taa = movea(t,ta); Paa = [taa._a,taa._b,taa._c]
            tab = moveb(t,ta); Pab = [tab._a,tab._b,tab._c]
            tac = movec(t,ta); Pac = [tac._a,tac._b,tac._c]

            tba = movea(t,tb); Pba = [tba._a,tba._b,tba._c]
            tbb = moveb(t,tb); Pbb = [tbb._a,tbb._b,tbb._c]
            tbc = movec(t,tb); Pbc = [tbc._a,tbc._b,tbc._c]

            tca = movea(t,tc); Pca = [tca._a,tca._b,tca._c]
            tcb = moveb(t,tc); Pcb = [tcb._a,tcb._b,tcb._c]
            tcc = movec(t,tc); Pcc = [tcc._a,tcc._b,tcc._c]

            P = vcat(P,Pa,Pb,Pc,Paa,Pab,Pac,Pba,Pbb,Pbc,Pca,Pcb,Pcc)
        end
        unique!(P)
    end

    inds = Array{Int}(undef,length(P))
    @inbounds for i ∈ eachindex(inds)
        temp = findfirst(isequal(P[i]),t.xy)
        inds[i] = isnothing(temp) ? 0 : temp
    end
    return inds[.!iszero.(inds)]
end


# extend VoronoiDelaunay methods to Tessellation
function VoronoiDelaunay.locate(t::Tessellation,x::Real,y::Real)
    @assert t.xmin ≤ x ≤ t.xmax && t.ymin ≤ y ≤ t.ymax "(x,y)=($x,$y) must be contained in ($(t.xmin),$(t.ymin)), ($(t.xmax),$(t.ymax))"
    return locate(t.tessellation,Point(tessellation_coords(x,y,t)...))
end
for fun ∈ (:movea, :moveb, :movec)
    @eval VoronoiDelaunay.$(fun)(tess::Tessellation,trig) = $(fun)(tess.tessellation,trig)
end


@recipe function f(t::Tessellation; delaunay=true, voronoi=false)
    aspect_ratio --> 1
    xlim --> (t.xmin,t.xmax)
    ylim --> (t.ymin,t.ymax)
    legend --> false

    if delaunay
        x, y = getplotxy(delaunayedges(t.tessellation))
        @series begin
            seriestype --> :path
            tessellation_inv_coords(x,y,t)
        end
    end
    if voronoi
        x, y = getplotxy(voronoiedges(t.tessellation))
        @series begin
            seriestype --> :path
            tessellation_inv_coords(x,y,t)
        end
    end
end


end

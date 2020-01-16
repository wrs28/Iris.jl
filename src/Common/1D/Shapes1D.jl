export Interval

using LinearAlgebra
using RecipesBase

"""
    Interval <: AbstractShape{1,2}

Fields: `start`, `stop`, `origin`, `a` (length of interval)
"""
struct Interval <: AbstractShape{1,2}
    start::Point{1,Cartesian}
    stop::Point{1,Cartesian}
    origin::Point{1,Cartesian}
    a::Float64
end

"""
    (::Interval)(::Point{1}) -> ::Bool

evaluates to true if the point is between `start` and `stop`.

Alternately, `(::Interval)(x,y,z...)`
"""
(i::Interval)(p::Point) = i.start ≤ p.x ≤ i.stop
(i::Interval)(x::Real,y...) = i(Point(x,y...))

"""
    Interval(start, stop; [reference=:center]) -> interval

1-Dimensional Interval.

`ref` determines where `origin` of the interval is relative to `start` & `stop`
"""
function Interval(start,stop;reference::Symbol=:center)
    p1, p2 = Point(start), Point(stop)
    return Interval(p1,p2,generate_origin(p1,p2,reference),norm(stop-start))
end

function generate_origin(start::Point{1},stop::Point{1},ref::Symbol)
    left,right = minmax(start,stop)
    if Base.sym_in(ref,(:left,:l,:Left,:L))
        origin = Point(left)
    elseif Base.sym_in(ref,(:right,:r,:Right,:R))
        origin = Point(right)
    else
        origin = Point((left+right)/2)
    end
    return origin
end


import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK

function Base.show(io::IO,i::Interval)
    printstyled(io, "Interval",color=PRINTED_COLOR_DARK)
    print(io," (")
    printstyled(io,i.start.x,color=PRINTED_COLOR_NUMBER)
    print(io,", ")
    printstyled(io,i.stop.x,color=PRINTED_COLOR_NUMBER)
    print(io,")")
end

################################################################################
# Plotting

@recipe function f(i::Interval)
    ylims --> (-1,1)
    legend --> false
    [i.start.x,i.stop.x],[0,0]
end

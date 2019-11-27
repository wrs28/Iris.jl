export Interval

using LinearAlgebra

"""
    struct Interval(start, stop; [ref=:center]) -> I

1-Dimensional Interval shape.

`(::Interval)(::Point)` evaluates to true if the point is between `start` and `stop`.
Alternately, `(::Interval)(x,y,z...)`

`ref` determines where `origin` of the interval is relative to `start` & `stop`
"""
struct Interval <: AbstractShape{1,2}
    start::Point{1}
    stop::Point{1}
    origin::Point{1}
    a::Float64
end
(i::Interval)(p::Point) = i.start ≤ p.x ≤ i.stop
(i::Interval)(x::Real,y...) = i(Point(x,y...))

function Interval(start,stop;ref::Symbol=:center)
    p1, p2 = Point(start), Point(stop)
    return Interval(p1,p2,generate_origin(p1,p2,ref),norm(stop-start))
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

function Base.show(io::IO,i::Interval)
    printstyled(io, "Interval",color=PRINTED_COLOR_DARK)
    print(io," (")
    printstyled(io,i.start,color=PRINTED_COLOR_NUMBER)
    print(io,",")
    printstyled(io,i.stop,color=PRINTED_COLOR_NUMBER)
    print(io,")")
end

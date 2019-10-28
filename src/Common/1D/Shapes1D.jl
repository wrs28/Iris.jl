export Interval


struct Interval <: AbstractShape{1,2}
    start::Float64
    stop::Float64
    origin::Point{1}
    a::Float64
end
(i::Interval)(p::Point) = i.start ≤ p.x ≤ i.stop
(i::Interval)(x::Real,y...) = i(Point(x,y...))

Interval(start,stop;ref::Symbol=:center) = Interval(start,stop,generate_origin(start,stop,ref),abs(stop-start))

function generate_origin(start::Real,stop::Real,ref::Symbol)
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

export Interval

using ..Points

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

# struct BraggStack{N,T}
#     origin::Float64
#     depths::NTuple{N,Float64}
#     indices::NTuple{N,T}
#     periods::Int
#     integrated_depths::NTuple{N,Float64}
#
#     function BraggStack(origin,depths,indices,periods)
#         new{length(depths),eltype(indices)}(origin,depths,indices,periods,tuple(cumsum(depths...)))
#     end
# end
# function BraggStack(p::Point)
#     0 ≤ (p.x-origin)
# end

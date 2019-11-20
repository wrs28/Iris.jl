"""
    module ElectricFields
"""
module ElectricFields

export ElectricField

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_INSTRUCTION
import ..PRINTED_COLOR_VARIABLE

using ..Points
using Interpolations

"""
    struct ElectricField{N}

    ElectricField(::Vector{Point},::Matrix) -> E
    ElectricField(::Vector{Point},::Int) -> E

Stores `M` electric fields. `E.x`, `E.y`, `E.z` create views into components of electric fields.

Access a given field (mode) `i` with `E(i)`.

The order in which the arguments are given to `ElectricField` don't matter.
"""
struct ElectricField{N} <:AbstractArray{ComplexF64,2}
    positions::Vector{Point{N}}
    values::Matrix{ComplexF64}
    starts::Vector{Float64}
    stops::Vector{Float64}
    start_inds::Vector{Int}
    stop_inds::Vector{Int}

    function ElectricField(pos::Vector{Point{N}}, val::AbstractMatrix) where N
        starts = Vector{Float64}(undef,N)
        stops = similar(starts)
        start_inds = Vector{Int}(undef,N)
        stop_inds = similar(start_inds)
        for i ∈ eachindex(starts)
            starts[i],start_inds[i] = findmin(map(p->p.vector[i],pos))
            stops[i],stop_inds[i] = findmax(map(p->p.vector[i],pos))
        end
        return new{N}(pos,val,starts,stops,start_inds,stop_inds)
    end
end
ElectricField(pos::Array{T},val::AbstractVector) where T<:Point = ElectricField(pos,reshape(val,length(val),1))
ElectricField(val,pos::Array{T}) where T<:Point = ElectricField(pos,val)
ElectricField(pos::Array{T,1},m::Int=1) where T<:Point = ElectricField(pos,zeros(ComplexF64,3length(pos),m))
ElectricField(m::Int,pos::Array{T,1}) where T<:Point = ElectricField(pos,m)

(e::ElectricField)(i) = ElectricField(e.pos,view(e.val,:,i))

Base.conj(e::ElectricField) = ElectricField(e.pos,conj(e.val))

fnames = (:+,:-)
for fname ∈ fnames
    @eval begin
        function Base.$fname(e1::ElectricField,e2::ElectricField)
            e1.pos == e2.pos || throw("cannot binary operate ElectricFields that do not share field points")
            return ElectricField(e1.pos,$fname(e1.val,e2.val))
        end
    end
end
Base.:*(a::Number,e::ElectricField) = ElectricField(e.pos,a*e.val)
Base.:\(a::Number,e::ElectricField) = ElectricField(e.pos,a\e.val)

Base.getindex(e::ElectricField,inds...) = e.val[inds...]
Base.setindex!(e::ElectricField,v,inds...) = setindex!(e.val,v,inds...)
fnames = (:firstindex, :lastindex, :length, :size, :eltype, :iterate)
for fname ∈ fnames @eval Base.$fname(e::ElectricField) = $fname(e.val) end
Base.iterate(e::ElectricField,state) = iterate(e.val,state)
Base.IndexStyle(::Type{<:ElectricField}) = IndexLinear()


function Base.getproperty(e::ElectricField{N},sym::Symbol) where N
    n,m = size(getfield(e,:values))
    if sym == :pos
        return getfield(e,:positions)
    elseif sym == :val
        return getfield(e,:values)
    elseif Base.sym_in(sym,(:x,:X))
        return view(getfield(e,:values),0n÷3+1:1n÷3,:)
    elseif Base.sym_in(sym,(:y,:Y))
        return view(getfield(e,:values),1n÷3+1:2n÷3,:)
    elseif Base.sym_in(sym,(:z,:Z))
        return view(getfield(e,:values),2n÷3+1:3n÷3,:)
    elseif N==1 && sym==:left
        start_ind = getfield(e,:start_inds)[1]
        values = getfield(e,:values)
        inds = [0n÷3+start_ind,1n÷3+start_ind,2n÷3+start_ind]
        return values[inds,:]
    elseif N==1 && sym==:right
        start_ind = getfield(e,:stop_inds)[1]
        values = getfield(e,:values)
        inds = [0n÷3+start_ind,1n÷3+start_ind,2n÷3+start_ind]
        return values[inds,:]
    elseif N==1 && sym==:start
        return getfield(e,:starts)[1]
    elseif N==1 && sym==:stop
        return getfield(e,:stops)[1]
    else
        return getfield(e,sym)
    end
end

function Base.propertynames(::ElectricField{1},private=false)
    if private
        return fieldnames(ElectricField)
    else
        return (:positions,:values,:x,:y,:z,:left,:right)
    end
end

# Pretty printing
Base.show(io::IO,mime::MIME"text/plain",e::ElectricField) = show(io,e)
function Base.show(io::IO,e::ElectricField)
    n,m = size(e)
    printstyled(io,"ElectricField",color=PRINTED_COLOR_DARK)
    print(io," (")
    if m > 1
        printstyled(io,m,color=PRINTED_COLOR_NUMBER)
        print(io," field")
        m > 1 ? print(io,"s @ ") : print(io," @ ")
    end
    printstyled(io,n÷3,color=PRINTED_COLOR_NUMBER)
    print(io," points)")
end


################################################################################
# extending INTERPOLATIONS

const AbInt = AbstractInterpolation
struct InterpolatedElectricField{N,TIX,TIY,TIZ}
    positions::Vector{Point{N}}
    x_interpolation::TIX
    y_interpolation::TIY
    z_interpolation::TIZ
    InterpolatedElectricField(pos::Vector{Point{N}},xi::TIX,yi::TIY,zi::TIZ) where {N,TIX<:AbInt,TIY<:AbInt,TIZ<:AbInt} =
            new{N,TIX,TIY,TIZ}(pos,xi,yi,zi)
end

function Base.getproperty(ie::InterpolatedElectricField,sym::Symbol)
    if sym==:x
        return getfield(ie,:x_interpolation)
    elseif sym==:y
        return getfield(ie,:y_interpolation)
    elseif sym==:z
        return getfield(ie,:z_interpolation)
    elseif sym==:pos
        return getfield(ie,:positions)
    else
        return getfield(ie,sym)
    end
end

Base.propertynames(::InterpolatedElectricField{N},private=false) where N = propertynames(ElectricField{N},private)

Interpolations.interpolate(e::ElectricField{1}, interpmode, gridstyle) =
    map(i->map(s->interpolate(getproperty(e(i),s), interpmode, gridstyle),(:x,:y,:z)),1:size(e,2))

fnames = (:LinearInterpolation,:CubicSplineInterpolation)
for fn ∈ fnames
    @eval begin
        function Interpolations.$fn(e::ElectricField{1}; kwargs...)
            x = map(p->p.x,e.pos)
            perm = sortperm(x)
            xs = LinRange(e.starts[1],e.stops[1],length(e.pos))
            itps = map(i->InterpolatedElectricField(e.pos,map(s->$fn(xs,getproperty(e(i),s)[perm]; kwargs...), (:x,:y,:z))...), ntuple(identity,size(e,2)))
            if length(itps)==1
                return itps[1]
            else
                return itps
            end
        end
    end
end

# Pretty Printing
function Base.show(io::IO,ie::InterpolatedElectricField{N}) where N
    print(io,"$(N)D ")
    printstyled(io,"Interpolated ElectricField ",color=PRINTED_COLOR_DARK)
    printstyled(io,"(access fields ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,".x",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,".y",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,".z",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
end

end #module

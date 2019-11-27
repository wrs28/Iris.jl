"""
Utilities for storing and manipulating electric fields
"""
module ElectricFields

export ElectricField

dimensional_files = (
    "1D/ElectricFields1D.jl",
    )

using ..Points
using Interpolations

"""
    ElectricField{N} <: AbstractMatrix{ComplexF64}

`E.x`, `E.y`, `E.z` create views into components of electric fields.

Access a given field (mode) `i` with `E(i)`.
"""
struct ElectricField{N} <: AbstractMatrix{ComplexF64}
    positions::Vector{Point{N}}
    values::Matrix{ComplexF64}
    start::Point{N}
    stop::Point{N}
    start_inds::Vector{Int}
    stop_inds::Vector{Int}

    function ElectricField(pos::Vector{Point{N}}, val::AbstractMatrix, start::Point{N}, stop::Point{N}, start_inds::Vector, stop_inds::Vector) where N
        size(val,1)==3length(pos) || throw("provided matrix has size(matrix,1)=$(size(val,1))≠3length(pos)=$(3length(pos))")
        return new{N}(pos,val,start,stop,start_inds,stop_inds)
    end
end

"""
    ElectricField(::Vector{Point},vals::Matrix) -> E

Store `vals` as electric field (`size(vals,1)` = number of points, `size(vals,2)` = number of fields)
"""
ElectricField(pos::Vector{T},val::AbstractVector,args...) where T<:Point = ElectricField(pos,reshape(val,length(val),1),args...)

"""
    ElectricField(::Vector{Point}, m::Integer) -> E

Stores `m` blank electric fields.
"""
ElectricField(arg,m::Integer) where T<:Point = ElectricField(arg,zeros(ComplexF64,3length(arg),m))

foreach(include,dimensional_files)

(e::ElectricField)(i) = ElectricField(e.pos,view(e.values,:,i),e.start,e.stop,e.start_inds,e.stop_inds)

Base.conj(e::ElectricField) = ElectricField(e.pos,conj(e.values))

fnames = (:+,:-)
for fname ∈ fnames
    @eval begin
        function Base.$fname(e1::ElectricField,e2::ElectricField)
            e1.pos == e2.pos || throw("cannot binary operate ElectricFields that do not share field points")
            return ElectricField(e1.pos,$fname(e1.val,e2.values))
        end
    end
end
Base.:*(a::Number,e::ElectricField) = ElectricField(e.pos,a*e.values)
Base.:*(e::ElectricField,a::Number) = a*e
Base.:\(a::Number,e::ElectricField) = ElectricField(e.pos,a\e.values)
Base.:/(e::ElectricField,a::Number) = a\e

Base.getindex(e::ElectricField,inds...) = e.values[inds...]
Base.setindex!(e::ElectricField,v,inds...) = setindex!(e.values,v,inds...)
fnames = (:firstindex, :lastindex, :length, :size, :eltype, :iterate)
for fname ∈ fnames @eval Base.$fname(e::ElectricField) = $fname(e.values) end
Base.iterate(e::ElectricField,state) = iterate(e.values,state)
Base.IndexStyle(::Type{<:ElectricField}) = IndexLinear()


function Base.getproperty(e::ElectricField{N}, sym::Symbol) where N
    p = getelectricfieldproperty(e,sym)
    if isnothing(p)
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
        else
            return getfield(e,sym)
        end
    else
        return p
    end
end

function Base.propertynames(::ElectricField{1},private=false)
    if private
        return fieldnames(ElectricField)
    else
        return (:positions,:values,:x,:y,:z,:left,:right)
    end
end


################################################################################
# extending INTERPOLATIONS

const AbInt = AbstractInterpolation

"""
Interpolated electric field.

Access interpolations through `E.x`, `E.y`, `E.z`. Evaluate with `E.x(::Real)`, etc.
"""
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

function Base.propertynames(::InterpolatedElectricField{N},private=false) where N
    if private
        return fieldnames(InterpolatedElectricField{N})
    else
        return propertynames(ElectricField{N})
    end
end

################################################################################
# Pretty Printing

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_INSTRUCTION
import ..PRINTED_COLOR_VARIABLE

function Base.show(io::IO,mime::MIME"text/plain",e::ElectricField)
    ds = displaysize(io)
    ioc = IOContext(io, :displaysize => (ds[1]÷3,ds[2]))
    show(io,e)
    printstyled(io,"\nX-component\n",color=PRINTED_COLOR_VARIABLE)
    show(ioc,mime,e.x)
    printstyled(io,"\nY-component\n",color=PRINTED_COLOR_VARIABLE)
    show(ioc,mime,e.y)
    printstyled(io,"\nZ-component\n",color=PRINTED_COLOR_VARIABLE)
    show(ioc,mime,e.z)
end

function Base.show(io::IO,e::ElectricField)
    n,m = size(e)
    printstyled(io,"ElectricField",color=PRINTED_COLOR_DARK)
    print(io," (")
    if m > 1
        printstyled(io,m,color=PRINTED_COLOR_NUMBER)
        print(io," field")
        m > 1 ? print(io,"s, ") : nothing
    end
    printstyled(io,n÷3,color=PRINTED_COLOR_NUMBER)
    print(io," points)")
end

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

"""
Utilities for storing and manipulating electric (vector) and scalar fields
"""
module VectorFields

export VectorField
export ScalarField
export ElectricField
export update!

dimensional_files = (
    "1D/VectorFields1D.jl",
    # "2D/VectorFields2D.jl",
    # "3D/VectorFields3D.jl",
    )

using ..Points
using Interpolations

"""
    VectorFields{N,M} <: AbstractMatrix{ComplexF64}

`N`is the dimension, `M` is the number of field components
"""
struct VectorField{N,M} <: AbstractMatrix{ComplexF64} # N = Dimension, M = number of field components, C = coordinate type
    positions::Vector{Point{N,Cartesian}}
    values::Matrix{ComplexF64}
    start::Point{N,Cartesian}
    stop::Point{N,Cartesian}
    start_inds::Vector{Int}
    stop_inds::Vector{Int}

    function VectorField{M}(pos::Vector{Point{N}}, val::AbstractMatrix, start::Point{N}, stop::Point{N}, start_inds::Vector, stop_inds::Vector) where {N,M}
        size(val,1)==M*length(pos) || throw("provided matrix has size(matrix,1)=$(size(val,1))≠$M*length(pos)=$(M*length(pos))")
        return new{N,M}(Cartesian.(pos),val,Cartesian(start),Cartesian(stop),start_inds,stop_inds)
    end
end

VectorField{M}(pos::Vector{TP},val::AbstractVector,args...) where {M,TP<:Point{N}} where N = VectorField{M}(pos,reshape(val,length(val),1),args...)

VectorField(f::VectorField{M},vals::AbstractVecOrMat) where M = VectorField{M}(f.pos,vals,f.start,f.stop,f.start_inds,f.stop_inds)

update!(f::VectorField,vals::AbstractVecOrMat) = copyto!(f.values,vals)

(f::VectorField)(i) = VectorField(f,view(f.values,:,i))

Base.conj(f::VectorField) = VectorField(f,conj(f.values))
Base.conj!(f::VectorField) = begin conj!(f.values); f end

fnames = (:+,:-)
for fname ∈ fnames
    @eval begin
        function Base.$fname(f1::VectorField{M},f2::VectorField{M}) where M
            f1.pos == f2.pos || throw("cannot binary operate `VectorField`s that do not share field points")
            return VectorField(f1,$fname(f1.values,f2.values))
        end
    end
end
Base.:*(f::VectorField,a::Number) = a*f
Base.:*(a::Number,f::VectorField) = VectorField(f,a*f.values)
Base.:/(f::VectorField,a::Number) = a\f
Base.:\(a::Number,f::VectorField) = VectorField(f,a\f.values)

Base.getindex(f::VectorField,inds...) = f.values[inds...]
Base.setindex!(f::VectorField,v,inds...) = setindex!(f.values,v,inds...)
fnames = (:firstindex, :lastindex, :length, :size, :eltype, :iterate)
for fname ∈ fnames @eval Base.$fname(f::VectorField) = $fname(f.values) end
Base.iterate(f::VectorField,state) = iterate(f.values,state)
Base.IndexStyle(::Type{<:VectorField}) = IndexLinear()

function Base.getproperty(e::VectorField{N,1}, sym::Symbol) where N
    p = getelectricfieldproperty(e,sym)
    if isnothing(p)
        n,m = size(getfield(e,:values))
        if sym == :pos
            return getfield(e,:positions)
        elseif sym == :val
            return getfield(e,:values)
        else
            return getfield(e,sym)
        end
    else
        return p
    end
end

function Base.getproperty(e::VectorField{N,3}, sym::Symbol) where N
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

"""
    ScalarField{N} <: AbstractMatrix{ComplexF64}

Access a given field (mode) `i` with `S(i)`.
"""
ScalarField{N} = VectorField{N,1}

"""
    ScalarField(pos, values, start, stop, start_inds, stop_inds) -> sf
"""
ScalarField(args...) = VectorField{1}(args...)

"""
    ScalarField(sf::ScalarField, values) -> sf
"""
ScalarField(f::ScalarField,vals::AbstractVecOrMat) = VectorField(f::ScalarField,vals::AbstractVecOrMat)

"""
    ElectricField{N} <: AbstractMatrix{ComplexF64}

`E.x`, `E.y`, `E.z` create views into components of electric fields.

Access a given field (mode) `i` with `E(i)`.
"""
ElectricField{N} = VectorField{N,3}

"""
    ElectricField{M}(pos, values, start, stop, start_inds, stop_inds) -> sf
"""
ElectricField(args...) = VectorField{3}(args...)

"""
    ElectricField(sf::ElectricField, values) -> sf
"""
ElectricField(f::ScalarField,vals::AbstractVecOrMat) = VectorField(f::ElectricField,vals::AbstractVecOrMat)


foreach(include,dimensional_files)

################################################################################
# extending INTERPOLATIONS
#
# const AbInt = AbstractInterpolation
#
# """
# Interpolated electric field.
#
# Access interpolations through `E.x`, `E.y`, `E.z`. Evaluate with `E.x(::Real)`, etc.
# """
# struct InterpolatedVectorField{N,M,TIX,TIY,TIZ}
#     positions::Vector{Point{N}}
#     x_interpolation::TIX
#     y_interpolation::TIY
#     z_interpolation::TIZ
#     InterpolatedVectorField{M}(pos::Vector{Point{N}},xi::TIX,yi::TIY,zi::TIZ) where {N,TIX<:AbInt,TIY<:AbInt,TIZ<:AbInt} =
#             new{N,M,TIX,TIY,TIZ}(pos,xi,yi,zi)
# end
#
# InterpolatedScalarField{N} = InterpolatedVectorField{N,1}
# InterpolatedVectorField{N} = InterpolatedVectorField{N,1}
#
# function Base.getproperty(ie::InterpolatedElectricField,sym::Symbol)
#     if sym==:x
#         return getfield(ie,:x_interpolation)
#     elseif sym==:y
#         return getfield(ie,:y_interpolation)
#     elseif sym==:z
#         return getfield(ie,:z_interpolation)
#     elseif sym==:pos
#         return getfield(ie,:positions)
#     else
#         return getfield(ie,sym)
#     end
# end
#
# function Base.propertynames(::InterpolatedElectricField{N},private=false) where N
#     if private
#         return fieldnames(InterpolatedElectricField{N})
#     else
#         return propertynames(ElectricField{N})
#     end
# end

################################################################################
# Pretty Printing

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_INSTRUCTION
import ..PRINTED_COLOR_VARIABLE

function Base.show(io::IO,mime::MIME"text/plain",f::ScalarField)
    ds = displaysize(io)
    ioc = IOContext(io, :displaysize => (ds[1],ds[2]))
    show(io,f)
    show(ioc,mime,v.x)
end

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

function Base.show(io::IO,f::VectorField{N,M}) where {N,M}
    n,m = size(f)
    if M==1
        printstyled(io,"Scalar",color=PRINTED_COLOR_DARK)
    elseif M==3
        printstyled(io,"ElectricField",color=PRINTED_COLOR_DARK)
    else
        printstyled(io,"$M-component VectorField",color=PRINTED_COLOR_DARK)
    end
    print(io," (")
    if m > 1
        printstyled(io,m,color=PRINTED_COLOR_NUMBER)
        print(io," field")
        m > 1 ? print(io,"s, ") : nothing
    end
    printstyled(io,n÷M,color=PRINTED_COLOR_NUMBER)
    print(io," points)")
end

# function Base.show(io::IO,ie::InterpolatedScalarField{N}) where N
#     print(io,"$(N)D ")
#     printstyled(io,"Interpolated ScalarField ",color=PRINTED_COLOR_DARK)
#     printstyled(io,"(access fields ",color=PRINTED_COLOR_INSTRUCTION)
#     printstyled(io,".values",color=PRINTED_COLOR_VARIABLE)
#     printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
# end
#
# function Base.show(io::IO,ie::InterpolatedElectricField{N}) where N
#     print(io,"$(N)D ")
#     printstyled(io,"Interpolated VectorField ",color=PRINTED_COLOR_DARK)
#     printstyled(io,"(access fields ",color=PRINTED_COLOR_INSTRUCTION)
#     printstyled(io,".x",color=PRINTED_COLOR_VARIABLE)
#     printstyled(io,", ",color=PRINTED_COLOR_INSTRUCTION)
#     printstyled(io,".y",color=PRINTED_COLOR_VARIABLE)
#     printstyled(io,", ",color=PRINTED_COLOR_INSTRUCTION)
#     printstyled(io,".z",color=PRINTED_COLOR_VARIABLE)
#     printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
# end

end #module

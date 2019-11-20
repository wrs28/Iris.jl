"""
    module InversionFields
"""
module InversionFields

export InversionField

using ...Common

import ...Common.PRINTED_COLOR_NUMBER
import ...Common.PRINTED_COLOR_DARK

"""
    struct InversionField{N}

    InversionField(::Vector{Point},::Vector) -> D

Stores an inversion field.

The order in which the arguments are given to `InversionField` don't matter.
"""
struct InversionField{N} <:AbstractVector{Float64}
    positions::Vector{Point{N}}
    values::Vector{Float64}
    InversionField(pos::Array{Point{N},1},val::AbstractVector) where N = new{N}(pos,val)
end
InversionField(val,pos::Array{T}) where T<:Point = InversionField(pos,val)
InversionField(pos::Array{T,1}) where T<:Point = InversionField(pos,zeros(Float64,length(pos)))

Base.conj(d::InversionField) = InversionField(d.pos,conj(d.val))

fnames = (:+,:-)
for fname ∈ fnames
    @eval begin
        function Base.$fname(d1::InversionField,d2::InversionField)
            d1.pos == d2.pos || throw("cannot binary operate InversionFields that do not share field points")
            return InversionField(d1.pos,$fname(d1.val,d2.val))
        end
    end
end
Base.:*(a::Number,d::InversionField) = InversionField(d.pos,a*d.val)
Base.:\(a::Number,d::InversionField) = InversionField(d.pos,a\d.val)

Base.getindex(d::InversionField,inds...) = d.val[inds...]
Base.getindex(d::InversionField,ind) = d.val[mod1(ind,length(d))]
Base.setindex!(d::InversionField,v,inds...) = setindex!(d.val,v,inds...)
fnames = (:firstindex, :lastindex, :length, :size, :eltype, :iterate)
for fname ∈ fnames @eval Base.$fname(d::InversionField) = $fname(d.val) end
Base.iterate(d::InversionField,state) = iterate(d.val,state)
Base.IndexStyle(::Type{<:InversionField}) = IndexLinear()


function Base.getproperty(d::InversionField,sym::Symbol)
    n = length(getfield(d,:values))
    if sym == :pos
        return getfield(d,:positions)
    elseif sym == :val
        return getfield(d,:values)
        else
        return getfield(d,sym)
    end
end

# Pretty printing
Base.show(io::IO,mime::MIME"text/plain",d::InversionField) = show(io,d)
function Base.show(io::IO,d::InversionField)
    n = length(d)
    printstyled(io,"InversionField",color=PRINTED_COLOR_DARK)
    print(io," @ ")
    printstyled(io,n,color=PRINTED_COLOR_NUMBER)
    print(io," points")
end

end #module

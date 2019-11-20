"""
    module PolarizationFields
"""
module PolarizationFields

export PolarizationField

using ...Common

import ...Common.PRINTED_COLOR_NUMBER
import ...Common.PRINTED_COLOR_DARK

"""
    struct PolarizationField{N}

    PolarizationField(::Vector{Point},::Vector) -> P

Stores a polarization field. `P.x`, `P.y`, `P.z` create views into components of polarization fields.

The order in which the arguments are given to `PolarizationField` don't matter.
"""
struct PolarizationField{N} <:AbstractVector{ComplexF64}
    positions::Vector{Point{N}}
    values::Vector{ComplexF64}
    PolarizationField(pos::Array{Point{N},1},val::AbstractVector) where N = new{N}(pos,val)
end
PolarizationField(val,pos::Array{T}) where T<:Point = PolarizationField(pos,val)
PolarizationField(pos::Array{T,1}) where T<:Point = PolarizationField(pos,zeros(ComplexF64,3length(pos)))

Base.conj(p::PolarizationField) = PolarizationField(p.pos,conj(p.val))

fnames = (:+,:-)
for fname ∈ fnames
    @eval begin
        function Base.$fname(p1::PolarizationField,p2::PolarizationField)
            p1.pos == p2.pos || throw("cannot binary operate PolarizationFields that do not share field points")
            return PolarizationField(p1.pos,$fname(p1.val,p2.val))
        end
    end
end
Base.:*(a::Number,p::PolarizationField) = PolarizationField(p.pos,a*p.val)
Base.:\(a::Number,p::PolarizationField) = PolarizationField(p.pos,a\p.val)

Base.getindex(p::PolarizationField,inds...) = p.val[inds...]
Base.setindex!(p::PolarizationField,v,inds...) = setindex!(p.val,v,inds...)
fnames = (:firstindex, :lastindex, :length, :size, :eltype, :iterate)
for fname ∈ fnames @eval Base.$fname(p::PolarizationField) = $fname(p.val) end
Base.iterate(p::PolarizationField,state) = iterate(p.val,state)
Base.IndexStyle(::Type{<:PolarizationField}) = IndexLinear()


function Base.getproperty(p::PolarizationField,sym::Symbol)
    n,m = size(getfield(p,:values))
    if sym == :pos
        return getfield(p,:positions)
    elseif sym == :val
        return getfield(p,:values)
    elseif Base.sym_in(sym,(:x,:X))
        return view(getfield(p,:values),0n÷3+1:1n÷3)
    elseif Base.sym_in(sym,(:y,:Y))
        return view(getfield(p,:values),1n÷3+1:2n÷3)
    elseif Base.sym_in(sym,(:z,:Z))
        return view(getfield(p,:values),2n÷3+1:3n÷3)
    else
        return getfield(p,sym)
    end
end

# Pretty printing
Base.show(io::IO,mime::MIME"text/plain",p::PolarizationField) = show(io,p)
function Base.show(io::IO,p::PolarizationField)
    n = length(p)
    printstyled(io,"PolarizationField",color=PRINTED_COLOR_DARK)
    print(io," @ ")
    printstyled(io,n÷3,color=PRINTED_COLOR_NUMBER)
    print(io," points")
end

end #module

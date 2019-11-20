"""
    module RealElectricFields
"""
module RealElectricFields

export RealElectricField

using ..Common
using RecipesBase

import ..Common.PRINTED_COLOR_NUMBER
import ..Common.PRINTED_COLOR_DARK

"""
    struct RealElectricField{N}

    RealElectricField(::Vector{Point},::Vector) -> D

Stores an inversion field.

The order in which the arguments are given to `InversionField` don't matter.
"""
struct RealElectricField{N} <:AbstractVector{Float64}
    positions::Vector{Point{N}}
    values::Vector{Float64}
    RealElectricField(pos::Array{Point{N},1},val::AbstractVector) where N = new{N}(pos,val)
end
RealElectricField(val,pos::Array{T}) where T<:Point = RealElectricField(pos,val)
RealElectricField(pos::Array{T,1}) where T<:Point = RealElectricField(pos,zeros(Float64,3length(pos)))

Base.conj(d::RealElectricField) = RealElectricField(d.pos,conj(d.val))

fnames = (:+,:-)
for fname ∈ fnames
    @eval begin
        function Base.$fname(e1::RealElectricField,e2::RealElectricField)
            e1.pos == e2.pos || throw("cannot binary operate RealElectricField that do not share field points")
            return RealElectricField(e1.pos,$fname(e1.val,e2.val))
        end
    end
end
Base.:*(a::Number,e::RealElectricField) = RealElectricField(e.pos,a*e.val)
Base.:\(a::Number,e::RealElectricField) = RealElectricField(e.pos,a\e.val)

Base.getindex(e::RealElectricField,inds...) = e.val[inds...]
Base.getindex(e::RealElectricField,ind) = e.val[mod1(ind,length(e))]
Base.setindex!(e::RealElectricField,v,inds...) = setindex!(e.val,v,inds...)
fnames = (:firstindex, :lastindex, :length, :size, :eltype, :iterate)
for fname ∈ fnames @eval Base.$fname(e::RealElectricField) = $fname(e.val) end
Base.iterate(e::RealElectricField,state) = iterate(e.val,state)
Base.IndexStyle(::Type{<:RealElectricField}) = IndexLinear()


function Base.getproperty(e::RealElectricField,sym::Symbol)
    n = length(getfield(e,:values))
    if sym == :pos
        return getfield(e,:positions)
    elseif sym == :val
        return getfield(e,:values)
    elseif Base.sym_in(sym,(:x,:X))
        return view(getfield(e,:values),0n÷3+1:1n÷3)
    elseif Base.sym_in(sym,(:y,:Y))
        return view(getfield(e,:values),1n÷3+1:2n÷3)
    elseif Base.sym_in(sym,(:z,:Z))
        return view(getfield(e,:values),2n÷3+1:3n÷3)
    else
        return getfield(d,sym)
    end
end

# Pretty printing
Base.show(io::IO,mime::MIME"text/plain",e::RealElectricField) = show(io,e)
function Base.show(io::IO,e::RealElectricField)
    n = length(e)
    printstyled(io,"RealElectricField",color=PRINTED_COLOR_DARK)
    print(io," @ ")
    printstyled(io,n÷3,color=PRINTED_COLOR_NUMBER)
    print(io," points")
end

@recipe function f(e::RealElectricField{1})
    ElectricField(e.pos,e.val)
end

end #module

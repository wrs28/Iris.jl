"""
    module ElectricFields
"""
module ElectricFields

export ElectricField

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK
using ..Points

struct ElectricField{N} <:AbstractArray{ComplexF64,2}
    positions::Vector{Point{N}}
    values::Matrix{ComplexF64}
    ElectricField(pos::Array{Point{N},1},val::AbstractMatrix) where N = new{N}(pos,val)
end
ElectricField(pos::Array{T},val::AbstractVector) where T<:Point = ElectricField(pos,reshape(val,length(val),1))
ElectricField(val,pos::Array{T}) where T<:Point = ElectricField(pos,val)
ElectricField(pos::Array{T,1},m::Int=1) where T<:Point = ElectricField(pos,Matrix{ComplexF64}(undef,3length(pos),m))
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


function Base.getproperty(e::ElectricField,sym::Symbol)
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
        m > 1 ? print(io,"s, ") : print(io,", ")
    end
    printstyled(io,n÷3,color=PRINTED_COLOR_NUMBER)
    print(io," points)")
end

end #module

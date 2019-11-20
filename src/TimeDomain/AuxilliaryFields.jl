"""
    module AuxilliaryFields
"""
module AuxilliaryFields

export AuxilliaryField

using ...Common

import Common.PRINTED_COLOR_NUMBER
import Common.PRINTED_COLOR_DARK

"""
    struct AuxilliaryField{N}

    AuxilliaryField(::Vector{Point},::Matrix) -> G
    AuxilliaryField(::Vector{Point},::Int) -> G

Stores `M` auxilliary fields. `G.x`, `G.y`, `G.z` create views into components of auxilliary fields.

Access a given field (mode) `i` with `G(i)`.

The order in which the arguments are given to `AuxilliaryField` don't matter.
"""
struct AuxilliaryField{N} <:AbstractArray{ComplexF64,2}
    positions::Vector{Point{N}}
    values::Matrix{ComplexF64}
    AuxilliaryField(pos::Array{Point{N},1},val::AbstractMatrix) where N = new{N}(pos,val)
end
AuxilliaryField(pos::Array{T},val::AbstractVector) where T<:Point = AuxilliaryField(pos,reshape(val,length(val),1))
AuxilliaryField(val,pos::Array{T}) where T<:Point = AuxilliaryField(pos,val)
AuxilliaryField(pos::Array{T,1},m::Int=1) where T<:Point = AuxilliaryField(pos,zeros(ComplexF64,3length(pos),m))
AuxilliaryField(m::Int,pos::Array{T,1}) where T<:Point = AuxilliaryField(pos,m)

(g::AuxilliaryField)(i) = AuxilliaryField(g.pos,view(g.val,:,i))

Base.conj(g::AuxilliaryField) = AuxilliaryField(g.pos,conj(g.val))

fnames = (:+,:-)
for fname ∈ fnames
    @eval begin
        function Base.$fname(g1::AuxilliaryField,g2::AuxilliaryField)
            p1.pos == g2.pos || throw("cannot binary operate ElectricFields that do not share field points")
            return AuxilliaryField(g1.pos,$fname(g1.val,g2.val))
        end
    end
end
Base.:*(a::Number,g::AuxilliaryField) = AuxilliaryField(g.pos,a*g.val)
Base.:\(a::Number,g::AuxilliaryField) = AuxilliaryField(g.pos,a\g.val)

Base.getindex(g::AuxilliaryField,inds...) = g.val[inds...]
Base.setindex!(g::AuxilliaryField,v,inds...) = setindex!(g.val,v,inds...)
fnames = (:firstindex, :lastindex, :length, :size, :eltype, :iterate)
for fname ∈ fnames @eval Base.$fname(g::AuxilliaryField) = $fname(g.val) end
Base.iterate(g::AuxilliaryField,state) = iterate(g.val,state)
Base.IndexStyle(::Type{<:AuxilliaryField}) = IndexLinear()


function Base.getproperty(g::AuxilliaryField,sym::Symbol)
    n,m = size(getfield(g,:values))
    if sym == :pos
        return getfield(g,:positions)
    elseif sym == :val
        return getfield(g,:values)
    elseif Base.sym_in(sym,(:x,:X))
        return view(getfield(g,:values),0n÷3+1:1n÷3,:)
    elseif Base.sym_in(sym,(:y,:Y))
        return view(getfield(g,:values),1n÷3+1:2n÷3,:)
    elseif Base.sym_in(sym,(:z,:Z))
        return view(getfield(g,:values),2n÷3+1:3n÷3,:)
    else
        return getfield(g,sym)
    end
end

# Pretty printing
Base.show(io::IO,mime::MIME"text/plain",g::AuxilliaryField) = show(io,g)
function Base.show(io::IO,g::AuxilliaryField)
    n,m = size(g)
    printstyled(io,"AuxilliaryField",color=PRINTED_COLOR_DARK)
    print(io," (")
    if m > 1
        printstyled(io,m,color=PRINTED_COLOR_NUMBER)
        print(io," field")
        m > 1 ? print(io,"s @ ") : print(io," @ ")
    end
    printstyled(io,n÷3,color=PRINTED_COLOR_NUMBER)
    print(io," points)")
end

end #module

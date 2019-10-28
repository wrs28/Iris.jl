module Curlcurls

export Curlcurl

files = (
    "1D/Curlcurls1D.jl",
    # "2D/Curlcurls2D.jl",
    # "3D/Curlcurls3D.jl"
    )

import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_VARIABLE
using ..Lattices
using LinearAlgebra
using SparseArrays

struct Curlcurl{N}
    cc0::SparseMatrixCSC{ComplexF64,Int}
    cc1::SparseMatrixCSC{ComplexF64,Int}
    cc2::SparseMatrixCSC{ComplexF64,Int}
end

foreach(include,files)

Base.conj(cc::Curlcurl{N}) where N = Curlcurl{N}(cc.cc0,-cc.cc1,cc.cc2)

function Base.show(io::IO,cc::Curlcurl{N}) where N
    print(io,"$(N)D ")
    printstyled(io,"Curlcurl",color=PRINTED_COLOR_DARK)
    print(io,"(call with ")
    printstyled(io,"ky",color=PRINTED_COLOR_VARIABLE)
    print(io,", ")
    printstyled(io,"kz",color=PRINTED_COLOR_VARIABLE)
    print(io,")")
end

end

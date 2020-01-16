"""
Creation and manipulation of ∇×∇× differential operator
"""
module Curlcurls

export Curlcurl

files = (
    "1D/Curlcurls1D.jl",
    "2D/Symmetric/Curlcurls2D.jl",
    # "3D/Curlcurls3D.jl"
    )


using ..Domains
using ..Lattices
using ..Laplacians
using ..Points
using SparseArrays

import LinearAlgebra: I
import ..Symmetric, ..Unsymmetric

struct Curlcurl{N,CLASS}
    cc0::SparseMatrixCSC{ComplexF64,Int} # constant
    cc1::SparseMatrixCSC{ComplexF64,Int} # linear in kx,ky,kz
    cc2::SparseMatrixCSC{ComplexF64,Int} # quadratic in kx,ky,kz
    complex_scaling::SparseMatrixCSC{ComplexF64,Int}
end

foreach(include,files)

function Base.conj!(cc::Curlcurl{N}) where N
    conj!(cc.cc0)
    conj!(cc.cc1)
    conj!(cc.cc2)
    return cc
end
Base.conj(cc::Curlcurl{N}) where N = Curlcurl{N}(conj(cc.cc0),conj(cc.cc1),conj(cc.cc2),conj(cc.complex_scaling))

################################################################################
# Pretty Printing

import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_VARIABLE

function Base.show(io::IO,cc::Curlcurl{N,CLASS}) where {N,CLASS}
    print(io,"$(N)D ")
    CLASS<:Symmetric ? print(io,"Symmetric ") : nothing
    CLASS<:Unsymmetric ? print(io,"Unsymmetric ") : nothing
    printstyled(io,"Curlcurl",color=PRINTED_COLOR_DARK)
    print(io,"(call with ")
    printstyled(io,"ky",color=PRINTED_COLOR_VARIABLE)
    print(io,", ")
    printstyled(io,"kz",color=PRINTED_COLOR_VARIABLE)
    print(io,")")
end

end

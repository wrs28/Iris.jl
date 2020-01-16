"""
Creation and manipulation of ∇² differential operator
"""
module Laplacians

export Laplacian

files = (
    "1D/Laplacians1D.jl",
    "2D/Symmetric/Laplacians2D.jl",
    # "3D/Laplacians3D.jl"
    )

using ..Lattices
using ..Points
using SparseArrays

import ..Symmetric, ..Unsymmetric
import LinearAlgebra: I

struct Laplacian{N,CLASS}
    l0::SparseMatrixCSC{ComplexF64,Int}
    complex_scaling::SparseMatrixCSC{ComplexF64,Int}
end

foreach(include,files)

function Base.conj!(lap::Laplacian{N}) where N
    conj!(lap.l0)
    return lap
end
Base.conj(lap::Laplacian{N,CLASS}) where {N,CLASS} = Laplacian{N,CLASS}(conj(lap.l0),lap.complex_scaling)

################################################################################
# PRETTY PRINTING

import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_VARIABLE

function Base.show(io::IO,::Laplacian{N,CLASS}) where {N,CLASS}
    print(io,"$(N)D ")
    CLASS<:Symmetric ? print(io,"Symmetric ") : nothing
    CLASS<:Unsymmetric ? print(io,"Unsymmetric ") : nothing
    printstyled(io,"Laplacian", color=PRINTED_COLOR_DARK)
end

end

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
using LinearAlgebra
using SparseArrays

struct Laplacian{N} l0::SparseMatrixCSC{ComplexF64,Int} end

foreach(include,files)

function Base.conj!(lap::Laplacian{N}) where N
    conj!(lap.l0)
    return lap
end
Base.conj(lap::Laplacian{N}) where N = Laplacian{N}(conj(lap.l0))

################################################################################
# PRETTY PRINTING

import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_VARIABLE

function Base.show(io::IO,::Laplacian{N}) where N
    print(io,"$(N)D ")
    printstyled(io,"Laplacian",color=PRINTED_COLOR_DARK)
end

end

"""
Utility for defining self-energies
"""
module SelfEnergies

export SelfEnergy

files = (
    "1D/SelfEnergies1D.jl",
    "2D/Symmetric/SelfEnergies2D.jl",
    # "3D/SelfEnergies3D.jl"
    )

using ..BoundaryConditions
using ..Domains
using ..Points
using SparseArrays

import ..Symmetric, ..Unsymmetric
import LinearAlgebra: I

struct SelfEnergy{N,CLASS,M,TF}
    Σ0::NTuple{M,SparseMatrixCSC{ComplexF64,Int}} # for terms independent of ky, kz (i.e. from 2nd derivatives)
    Σ1::NTuple{M,SparseMatrixCSC{ComplexF64,Int}} # for terms linear in ky, kz (i.e. from 1st derivatives)
    Σ2::SparseMatrixCSC{ComplexF64,Int} # for terms quadratic in ky, kz (i.e. no derivatives)
    f::TF
end

foreach(include,files)

################################################################################
# Pretty Printing
import ..PRINTED_COLOR_DARK

function Base.show(io::IO,Σ::SelfEnergy{N,CLASS,M}) where {N,CLASS,M}
    print(io,"$(N)D ")
    CLASS<:Symmetric ? print(io,"Symmetric ") : nothing
    CLASS<:Unsymmetric ? print(io,"Unsymmetric ") : nothing
    printstyled(io,"SelfEnergy",color=PRINTED_COLOR_DARK)
    print(io," with ", M," components")
end

end # module

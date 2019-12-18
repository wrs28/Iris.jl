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

struct SelfEnergy{N,CLASS,M,TF}
    Σ0::NTuple{M,SparseMatrixCSC{ComplexF64,Int}} # for terms independent of ky, kz
    Σ1::NTuple{M,SparseMatrixCSC{ComplexF64,Int}} # for terms linear in ky, kz
    Σ2::SparseMatrixCSC{ComplexF64,Int} # for terms quadratic in ky, kz
    f::TF
end

foreach(include,files)

################################################################################
# Pretty Printing
import ..PRINTED_COLOR_DARK

function Base.show(io::IO,Σ::SelfEnergy{N,CLASS,M}) where {N,CLASS,M}
    print(io,"$(N)D $CLASS ")
    printstyled(io,"SelfEnergy",color=PRINTED_COLOR_DARK)
    print(io," with ", M," components")
end

end # module

"""
Utility for defining self-energies
"""
module SelfEnergies

export SelfEnergy

files = (
    "1D/SelfEnergies1D.jl",
    # "2D/SelfEnergies2D.jl",
    # "3D/SelfEnergies3D.jl"
    )

using ..BoundaryConditions
using ..Domains
using SparseArrays

struct SelfEnergy{N,M,TF}
    Σ0::NTuple{M,SparseMatrixCSC{ComplexF64,Int}}
    Σ1::NTuple{M,SparseMatrixCSC{ComplexF64,Int}}
    Σ2::SparseMatrixCSC{ComplexF64,Int}
    f::TF
end

foreach(include,files)

################################################################################
# Pretty Printing
import ..PRINTED_COLOR_DARK

function Base.show(io::IO,Σ::SelfEnergy{N,M}) where {N,M}
    print(io,"$(N)D ")
    printstyled(io,"SelfEnergy",color=PRINTED_COLOR_DARK)
    print(io," with ", M," components")
end

end # module

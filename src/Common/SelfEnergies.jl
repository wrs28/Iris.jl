module SelfEnergies

export SelfEnergy

using ..BoundaryConditions
using ..Domains
using SparseArrays

struct SelfEnergy{N,M,TF}
    Σ0::NTuple{M,SparseMatrixCSC{ComplexF64,Int}}
    Σ1::NTuple{M,SparseMatrixCSC{ComplexF64,Int}}
    Σ2::SparseMatrixCSC{ComplexF64,Int}
    f::TF
end

include("1D/SelfEnergies.jl")
# include("2D/SelfEnergies.jl")
# include("3D/SelfEnergies.jl")

function Base.show(io::IO,Σ::SelfEnergy{N,M}) where {N,M}
    print(io,"SelfEnergy in $(N)D with $(M) parts")
end

end

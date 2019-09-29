module Curlcurls

export Curlcurl

using ..Lattices
using LinearAlgebra
using SparseArrays

struct Curlcurl{N}
    cc0::SparseMatrixCSC{ComplexF64,Int}
    cc1::SparseMatrixCSC{ComplexF64,Int}
    cc2::SparseMatrixCSC{ComplexF64,Int}
end

include("1D/Curlcurls.jl")
# include("2D/Curlcurls.jl")
# include("3D/Curlcurls.jl")

function Base.show(io::IO,cc::Curlcurl{N}) where N
    print(io,"Curlcurl in $(N)D")
end

end

"""
    module MaxwellOperators
"""
module MaxwellOperators

export Maxwell
export maxwell

files = (
    "1D/MaxwellOperators1D.jl",
    # "2D/MaxwellOperators2D.jl",
    # "3D/MaxwellOperators3D.jl"
    )

import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_VARIABLE

using ..ElectricFields
using ..Simulations
using ..Dispersions
using LinearAlgebra
using SparseArrays

const I3 = sparse(I,3,3)

struct Maxwell{N,TΣ,TSIM}
    A::SparseMatrixCSC{ComplexF64,Int}
    D²::SparseMatrixCSC{ComplexF64,Int}
    αεpχ::SparseMatrixCSC{Complex{Float64},Int}
    curlcurl::SparseMatrixCSC{ComplexF64,Int}
    Σs::TΣ
    αε::SparseMatrixCSC{ComplexF64,Int}
    αχ::SparseMatrixCSC{ComplexF64,Int}
    αχs::Vector{SparseMatrixCSC{ComplexF64,Int}}
    αdχdψr::Vector{SparseMatrixCSC{ComplexF64,Int}}
    αdχdψi::Vector{SparseMatrixCSC{ComplexF64,Int}}
    sim::TSIM
    kx::Float64
    ky::Float64
    kz::Float64
    ka::Float64
    kb::Float64
    kc::Float64
end


# load 1D, 2D, 3D
foreach(include,files)


"""
    maxwell(simulation{1}; ky=0, kz=0)
"""
maxwell(args...;kwargs...) = Maxwell(args...;kwargs...)


function Base.show(io::IO,::Maxwell{N}) where N
    print(io,N,"D ")
    printstyled(io,"Maxwell ",color=PRINTED_COLOR_DARK)
    print(io,"Operator (call w/ args ")
    printstyled(io,"ω",color=PRINTED_COLOR_VARIABLE)
    print(io,", [")
    printstyled(io,"ωs",color=PRINTED_COLOR_VARIABLE)
    print(io,",")
    printstyled(io,"ψs",color=PRINTED_COLOR_VARIABLE)
    print(io,"]) -> ")
    printstyled(io,"A",color=PRINTED_COLOR_VARIABLE)
end

end #module

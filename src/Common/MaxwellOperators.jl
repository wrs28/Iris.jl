"""
    module MaxwellOperators
"""
module MaxwellOperators

export Maxwell

files = (
    "1D/MaxwellOperators1D.jl",
    # "2D/MaxwellOperators2D.jl",
    # "3D/MaxwellOperators3D.jl"
    )

import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_VARIABLE
import ..PRINTED_COLOR_INSTRUCTION

using ..VectorFields
using ..Simulations
using ..Dispersions
using LinearAlgebra
using SparseArrays

const I3 = sparse(I,3,3)

"""
    maxwell(simulation{1}; ky=0, kz=0) -> m
"""
struct Maxwell{N,TΣ,TSIM}
    A::SparseMatrixCSC{ComplexF64,Int}
    D²::SparseMatrixCSC{ComplexF64,Int}
    αεpFχ::SparseMatrixCSC{Complex{Float64},Int}
    curlcurl::SparseMatrixCSC{ComplexF64,Int}
    Σs::TΣ
    αε::SparseMatrixCSC{ComplexF64,Int}
    Fχ::SparseMatrixCSC{ComplexF64,Int}
    Fχs::Vector{SparseMatrixCSC{ComplexF64,Int}}
    dFχdψr::Array{ComplexF64,3}
    dFχdψi::Array{ComplexF64,3}
    dFχdω::Array{ComplexF64,3}
    dFχdϕ::Array{ComplexF64,3}
    simulation::TSIM
    kx::Float64
    ky::Float64
    kz::Float64
    ka::Float64
    kb::Float64
    kc::Float64
end

# load 1D, 2D, 3D
foreach(include,files)


# Pretty Printing
function Base.show(io::IO,::Maxwell{N}) where N
    print(io,N,"D ")
    printstyled(io,"Maxwell ",color=PRINTED_COLOR_DARK)
    print(io,"Operator ")
    printstyled(io,"(call w/ args ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ω",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", [",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ωs",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,",",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ψs",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,"]) -> ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"A",color=PRINTED_COLOR_VARIABLE)
end

function Base.getproperty(m::Maxwell,sym::Symbol)
    if sym==:sim
        return getfield(m,:simulation)
    else
        return getfield(m,sym)
    end
end

function Base.propertynames(::Maxwell{1},private=false)
    if private
        return fieldnames(Maxwell)
    else
        (:simulation,:ky,:kz)
    end
end

end #module

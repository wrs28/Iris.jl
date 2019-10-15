module Scattering

export scattering
export green_fn

using ..Common
using ..Spectral
using ArnoldiMethodTransformations
using LinearAlgebra
using SparseArrays

include("1D/Scattering.jl")

function scattering_core(
            A::SparseMatrixCSC,
            j::AbstractArray)

    y = Array{ComplexF64}(undef,length(j))
    scattering_core!(y,A,j)
    return y
end
function scattering_core!(
            y::AbstractArray,
            A::SparseMatrixCSC,
            j::AbstractArray)

    length(j)==size(A,1) || throw("length(j)=$(length(j)) must equal size of system $(size(A,1))")
    A_lu! = ArnoldiMethodWrapper.ShiftAndInvert(A,0)
    A_lu!(y,j)
    return nothing
end

end # module

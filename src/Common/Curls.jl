"""
    module Curls
"""
module Curls

export Curl

files = (
    "1D/Curls1D.jl",
    # "2D/Curls2D.jl",
    # "3D/Curls3D.jl"
    )

using LinearAlgebra
using SparseArrays
using ..Lattices

import ..PRINTED_COLOR_DARK
import ..PRINTED_COLOR_VARIABLE
import ..PRINTED_COLOR_INSTRUCTION

struct Curl{N} curl::SparseMatrixCSC{Float64,Int} end

foreach(include,files)

# Pretty Printing
function Base.show(io::IO,cc::Curl{N}) where N
    print(io,"$(N)D ")
    printstyled(io,"Curl",color=PRINTED_COLOR_DARK)
    printstyled(io,"(call with no args or ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"ky",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,", ",color=PRINTED_COLOR_INSTRUCTION)
    printstyled(io,"kz",color=PRINTED_COLOR_VARIABLE)
    printstyled(io,")",color=PRINTED_COLOR_INSTRUCTION)
end

end

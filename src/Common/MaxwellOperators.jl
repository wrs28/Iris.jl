"""
    module MaxwellOperators
"""
module MaxwellOperators

export maxwell
export maxwell_lep
export maxwell_nep

using ..Simulations
using ..Dispersions
using Flatten
using LinearAlgebra
using NonlinearEigenproblems
using SparseArrays

const I3 = sparse(I,3,3)

struct Maxwell_LEP{N,TSIM,TΣ}
    sim::TSIM
    Σs::TΣ
    curlcurl::SparseMatrixCSC{ComplexF64,Int}
    αε::SparseMatrixCSC{Complex{Float64},Int}
    αχ::SparseMatrixCSC{Complex{Float64},Int}
    kx::Float64
    ky::Float64
    kz::Float64
    ka::Float64
    kb::Float64
    kc::Float64

    function Maxwell_LEP(
                sim::Simulation{N},
                Σs::Tuple,
                curlcurl::SparseMatrixCSC{ComplexF64,Int},
                αε::SparseMatrixCSC{ComplexF64,Int},
                αχ::SparseMatrixCSC{Complex{Float64},Int},
                kx::Number,
                ky::Number,
                kz::Number,
                ka::Number,
                kb::Number,
                kc::Number) where N

        return new{N,typeof(sim),typeof(Σs)}(sim,Σs,curlcurl,αε,αχ,kx,ky,kz,ka,kb,kc)
    end
end
function Base.show(io::IO,::Maxwell_LEP{N}) where N
    print(io,"Maxwell Linear Eigenproblem (call w/ args ")
    printstyled(io,"ω",color=:cyan)
    print(io,", [")
    printstyled(io,"ωs",color=:cyan)
    print(io,",")
    printstyled(io,"ψs",color=:cyan)
    print(io,"] -> ")
    printstyled(io,"A",color=:cyan)
    print(io,", ")
    printstyled(io,"B",color=:cyan)
    print(io,")")
end

struct Maxwell{N,TM}
    lep::TM
    Maxwell(lep::Maxwell_LEP{N}) where N = new{N,typeof(lep)}(lep)
end
function Base.show(io::IO,::Maxwell{N}) where N
    print(io,"Maxwell Operator (call w/ args ")
    printstyled(io,"ω",color=:cyan)
    print(io,", [")
    printstyled(io,"ωs",color=:cyan)
    print(io,",")
    printstyled(io,"ψs",color=:cyan)
    print(io,"])")
end

include("1D/MaxwellOperators.jl")
# include("2D/MaxwellOperators.jl")
# include("3D/MaxwellOperators.jl")

"""
    maxwell_lep(simulation; ky=0, kz=0)
"""
maxwell_lep

"""
    maxwell(simulation; ky=0, kz=0)
"""
maxwell

"""
    maxwell_nep(simulation; ky=0, kz=0, kwargs...)
"""
maxwell_nep

end #module

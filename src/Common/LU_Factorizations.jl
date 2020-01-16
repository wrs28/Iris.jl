"""
Utility for performing LU factorizations with UMFPACK, Pardiso, or MUMPS
"""
module LU_Factorizations

export LUfact
export MSolver
export PSolver
export USolver
export DEFAULT_LUPACK
export AbstractLUPACK

using ArnoldiMethodTransformations
using ArnoldiMethodTransformations: AbstractSolver
using LinearAlgebra
# using MPI
# using MUMPS
using Pardiso
using SparseArrays

import ..LUPACK
@eval const DEFAULT_LUPACK = $(LUPACK)()
"""
    DEFAULT_LUPACK = $DEFAULT_LUPACK
"""
DEFAULT_LUPACK

"""
    MSolver
"""
MSolver

"""
    PSolver
"""
PSolver

"""
    USolver
"""
USolver

"""
    AbstractLUPACK
"""
AbstractLUPACK = AbstractSolver

# colors for pretty printing
import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_LIGHT

# main structure used internally
struct LUfact{M,TLU,T}
    LU::TLU
    Z::SparseMatrixCSC{T,Int}
end

LinearAlgebra.lu(A,::T,args...;kwargs...) where T<:AbstractSolver  = lu(A,T,args...;kwargs...)
function LinearAlgebra.lu(A::AbstractArray{T},::Type{USolver},args...;kwargs...) where T
    alu = lu(A,args...;kwargs...)
    return LUfact{USolver,typeof(alu),T}(alu,spzeros(T,size(A)...))
end
function LinearAlgebra.lu(A::AbstractArray{T},::Type{MSolver},args...;kwargs...) where T
    throw("MUMPS still under construction")
    try MPI catch; throw(ErrorException("MUMPS not loaded. Try again after `using MUMPS`")) end
    MPI.Initialized() ? nothing : MPI.Init()
    alu = mumps_factorize(α)
    return LUfact{MSolver,typeof(alu),T}(alu,spzeros(T,size(A)...))
end
function LinearAlgebra.lu(A::AbstractArray{T},::Type{PSolver},args...;kwargs...) where T
    try PardisoSolver catch; throw(ErrorException("Pardiso not loaded. Try again after `using Pardiso`")) end
    issym = issymmetric(A)
    type = eltype(A)
    alu = PardisoSolver()
    set_iparm!(alu,1,1) # don't revert to defaults
    set_iparm!(alu,12,1) # transpose b/c of CSR vs SCS
    x = Vector{type}(undef,size(A,1))
    y = Vector{type}(undef,size(A,1))
    if issym & (type<:Real)
        set_matrixtype!(alu,Pardiso.REAL_SYM)
        pardiso(alu,x,sparse(tril(A)),y)
    elseif issym & (type<:Complex)
        set_matrixtype!(alu,Pardiso.COMPLEX_SYM)
        pardiso(alu,x,sparse(tril(A)),y)
    elseif !issym & (type<:Real)
        set_matrixtype!(alu,Pardiso.REAL_NONSYM)
        pardiso(alu,x,sparse(A),y)
    else # if !issym & type<:Complex
        set_matrixtype!(alu,Pardiso.COMPLEX_NONSYM)
        pardiso(alu,x,sparse(A),y)
    end
    set_phase!(alu,12) # analyze and factorize
    set_phase!(alu,33) # set to solve for future calls
    return LUfact{PSolver,typeof(alu),T}(alu,spzeros(T,size(A)...))
end

# define action of :LUfact
LinearAlgebra.ldiv!(Y,A::LUfact{USolver},B) = ldiv!(Y, A.LU, B)
LinearAlgebra.ldiv!(A::LUfact{USolver},B) = ldiv!(A.LU, B)
LinearAlgebra.ldiv!(Y,A::LUfact{MSolver},B) = ldiv!(Y, A.LU, B)
LinearAlgebra.ldiv!(Y,A::LUfact{PSolver},B) = pardiso(A.LU,Y,A.Z,B)
function Base.:\(A::LUfact,B)
    Y = similar(B)
    ldiv!(Y,A,B)
    return Y
end

# pretty printing
function Base.show(io::IO,si::LUfact{M}) where M
    printstyled(io,"LUfact ",color=PRINTED_COLOR_LIGHT)
    print(io,"(using ")
    if M<:MSolver
        printstyled(io,"MUMPS",color=PRINTED_COLOR_NUMBER)
    elseif M<:PSolver
        printstyled(io,"Pardiso",color=PRINTED_COLOR_NUMBER)
    elseif M<:USolver
        printstyled(io,"UMFPACK",color=PRINTED_COLOR_NUMBER)
    end
    print(io,")")
end

end # module

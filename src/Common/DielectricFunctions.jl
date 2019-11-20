#TODO: make ε a tensor
"""
    module DielectricFunctions

for defining dielectric and pump functions either as piecewise constant
or with user-supplied functions
"""
module DielectricFunctions

export DielectricFunction
export PumpFunction
export piecewise_constant_ε
export piecewise_constant_F

import ..PRINTED_COLOR_NUMBER
import ..PRINTED_COLOR_DARK

using ..Points

"""
    struct DielectricFunction

---

    DielectricFunction(n1=1,n2=0) -> dielectric
    DielectricFunction(n::Complex) -> dielectric

piecewise constant dielectric function with index given by `n1+1im*n2`

---

    DielectricFunction(ε::Function, parameters::Dict) -> dielectric

ε is any function with the signature ε(x,y,parameters::Dict{Symbol,Float64})
which evaluates to the dielectric at position `x,y`. For example:
    ε(x,y,params) = 2+sin(2πx/params[:period])

---

    (::DielectricFunction)(x,y) -> ε
"""
struct DielectricFunction{TF}
    ε::TF
    parameters::Dict{Symbol,Float64}

    function DielectricFunction(ε::TF, parameters::Dict{Symbol,T}) where {T<:Number,TF<:Function}
        if TF<:typeof(piecewise_constant_ε)
            if haskey(parameters,:n₁)
                get!(parameters,:n1,parameters[:n₁])
            elseif haskey(parameters,:n1)
                get!(parameters,:n₁,parameters[:n1])
            else
                throw(ErrorException("Dict parameters must contain at least one of :n₁, :n1"))
            end
            if haskey(parameters,:n₂)
                get!(parameters,:n2,parameters[:n₂])
            elseif haskey(parameters,:n2)
                get!(parameters,:n₂,parameters[:n2])
            else
                throw(ErrorException("Dict parameters must contain at least one of :n₂, :n2"))
            end
        end
        return new{TF}(ε,parameters)
    end
end
DielectricFunction() = DielectricFunction(piecewise_constant_ε,Dict(:n₁=>1,:n₂=>0))
DielectricFunction(n1::Real,n2::Real=0.0) = DielectricFunction(piecewise_constant_ε, Dict{Symbol,Float64}(:n₁=>n1,:n₂=>n2))
DielectricFunction(n::Complex) = DielectricFunction(real(n),imag(n))
DielectricFunction(ε::Function,n1::Real,n2::Real=0.0) = DielectricFunction(ε, Dict{Symbol,Float64}(:n₁=>n1,:n₂=>n2))
(df::DielectricFunction)(p::Point) = df.ε(p,df.parameters)

Base.getindex(de::DielectricFunction,sym::Symbol) = Base.getproperty(de,sym)
function Base.getproperty(de::DielectricFunction,sym::Symbol)
    if Base.sym_in(sym,(:n1,:n₁))
        return getfield(de,:parameters)[:n₁]
    elseif Base.sym_in(sym,(:n2,:n₂))
        return getfield(de,:parameters)[:n₂]
    else
        getfield(de,sym)
    end
end

Base.propertynames(::DielectricFunction,private=false) = (:n₁,:n₂,:parameters)

# Pretty Printing
function Base.show(io::IO,d::DielectricFunction)
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    printstyled(io, "DielectricFunction ",color=PRINTED_COLOR_DARK)
    if d.ε==piecewise_constant_ε
        print(io, "piecewise_constant_ε (n = ")
        printstyled(io,d.parameters[:n₁]+1im*d.parameters[:n₂],color=PRINTED_COLOR_NUMBER)
        print(io,")")
    else
        print(io,"(")
        count = 0
        for p ∈ d.parameters
            count += 1
            print(io,p.first,"=")
            printstyled(io,p.second,color=PRINTED_COLOR_NUMBER)
            count<length(d.parameters) ? print(io,", ") : nothing
        end
        print(io,")")
    end
end


"""
    struct PumpFunction

---

    PumpFunction(F=0) -> pump

piecewise constant pump function with pump parameter given by `F`

---

    PumpFunction(F::Function, parameters::Dict) -> pump

F is any real function with the signature F(x,y,parameters::Dict{Symbol,Float64})
which evaluates to the pump profile at position `x,y`. For example:
    F(x,y,params) = sin(2πx/params[:period])

---

    (::PumpFunction)(x,y) -> F
"""
struct PumpFunction{TF}
    F::TF
    parameters::Dict{Symbol,Float64}

    PumpFunction(F::TF, parameters::Dict) where TF<:Function = new{TF}(F,parameters)
end
PumpFunction() = PumpFunction(piecewise_constant_F,Dict{Symbol,Float64}(:F=>0))
PumpFunction(F::Real) = PumpFunction(piecewise_constant_F,Dict(:F=>F))
(pf::PumpFunction)(p::Point) = pf.F(p,pf.parameters)

# Pretty Printing
function Base.show(io::IO,d::PumpFunction)
    get(io,:tabbed2,false) ? print(io,"\t\t") : nothing
    printstyled(io, "PumpFunction ",color=PRINTED_COLOR_DARK)
    if d.F==piecewise_constant_F
        print(io, "piecewise_constant_F (F = ")
        printstyled(io, d.parameters[:F],color=PRINTED_COLOR_NUMBER)
        print(io, ")")
    else
        count = 0
        print(io,"(")
        for p ∈ d.parameters
            count += 1
            print(io,p.first,"=")
            printstyled(io, p.second,color=PRINTED_COLOR_NUMBER)
            count<length(d.parameters) ? print(io,", ") : nothing
        end
        print(io,")")
    end
end


"""
    piecewise_constant_ε(point,parameters) -> ε
"""
function piecewise_constant_ε(args...)
    parameters = args[end]
    n1 = get!(parameters,:n₁,nothing)
    n1 = get!(parameters,:n1,nothing)
    n2 = get!(parameters,:n₂,n1)
    n2 = get!(parameters,:n2,n2)
    return complex(n1,n2)^2
end


"""
    piecewise_constant_F(point,parameters) -> F
"""
piecewise_constant_F(args...) = get!(args[end],:F,nothing)

Base.conj(de::DielectricFunction{typeof(piecewise_constant_ε)}) = DielectricFunction(de.n1,-de.n2)

end # module

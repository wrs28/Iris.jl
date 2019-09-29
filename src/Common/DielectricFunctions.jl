#TODO: make ε a tensor
"""
    module DielectricFunctions

for defining dielectric and pump functions either as piecewise constant
or with user-supplied functions
"""
module DielectricFunctions

using ..Points

export DielectricFunction,
PumpFunction,
piecewise_constant_ε,
piecewise_constant_F


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
        return de.parameters[:n₁]
    elseif Base.sym_in(sym,(:n2,:n₂))
        return de.parameters[:n₂]
    else
        getfield(de,sym)
    end
end

function Base.show(io::IO,d::DielectricFunction)
    println(io, "DielectricFunction: ")
    if d.ε==piecewise_constant_ε
        print(io, "\tpiecewise_constant_ε: n = ",d.parameters[:n₁]+1im*d.parameters[:n₂])
    else
        for p ∈ d.parameters
            print(io, "\t", p.first, ": ")
            println(io, p.second)
        end
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

function Base.show(io::IO,d::PumpFunction)
    println(io, "PumpFunction: ")
    if d.F==piecewise_constant_F
        print(io, "\tpiecewise_constant_F: ",d.parameters[:F])
    else
        for p ∈ d.parameters
            print(io, "\t", p.first, ": ")
            println(io, p.second)
        end
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


end # module

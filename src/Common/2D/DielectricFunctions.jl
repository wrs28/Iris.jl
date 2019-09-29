"""
    module DielectricFunctions

for defining dielectric and pump functions either as piecewise constant
or with user-supplied functions
"""
module DielectricFunctions

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
struct DielectricFunction
    ε::Function
    parameters::Dict{Symbol,Float64}

    DielectricFunction(ε::Function=piecewise_constant_ε,
        parameters::Dict{Symbol,T}=Dict(:n₁=>1,:n₂=>0)
        ) where T<:Number = new(ε,parameters)

    DielectricFunction(n::Complex) = DielectricFunction(real(n),imag(n))

    function DielectricFunction(n1::Real,n2::Real=0.0)
        parameters = Dict(:n₁=>n1,:n₂=>n2)
        return DielectricFunction(piecewise_constant_ε,parameters)
    end

    (df::DielectricFunction)(x::Number,y::Number) = df.ε(x,y,df.parameters)

    function Base.show(io::IO,d::DielectricFunction)
        print(io, "DielectricFunction: ")
        if d.ε==piecewise_constant_ε
            print(io, "piecewise_constant_ε: ",d.parameters[:n₁]+1im*d.parameters[:n₂])
        else
            print(io, d.ε,": ", d.parameters)
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
struct PumpFunction
    F::Function
    parameters::Dict{Symbol,Float64}

    PumpFunction(F::Function=piecewise_constant_F,
        parameters::Dict{Symbol,T}=Dict(:F=>0)
        ) where T<:Number = new(F,parameters)

    function PumpFunction(F::Real)
        parameters = Dict(:F=>F)
        return PumpFunction(piecewise_constant_F,parameters)
    end

    (pf::PumpFunction)(x::Number,y::Number) = pf.F(x,y,pf.parameters)

    function Base.show(io::IO,d::PumpFunction)
        print(io, "PumpFunction: ")
        if d.F==piecewise_constant_F
            print(io, "piecewise_constant_F: ",d.parameters[:F])
        else
            print(io, d.F,": ", d.parameters)
        end
    end
end


"""
    piecewise_constant_ε(x,y,parameters) -> ε
"""
function piecewise_constant_ε(x::Number,y::Number,parameters::Dict)
    n1 = get!(parameters,:n₁,nothing)
    n2 = get!(parameters,:n₂,nothing)
    return complex(n1,n2)^2
end


"""
    piecewise_constant_F(x,y,parameters) -> F
"""
piecewise_constant_F(x::Number,y::Number,parameters::Dict) = get!(parameters,:F,nothing)


end # module

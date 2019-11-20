"""
    module MaxwellFDTD
"""
module MaxwellFDTD

export MaxwellBloch
export MaxwellBlochField
export maxwell_bloch

using LinearAlgebra
using ..Common
using ..PolarizationFields
using ..InversionFields
using ..AuxilliaryFields

import ..Common.PRINTED_COLOR_DARK
import ..Common.PRINTED_COLOR_NUMBER

struct MaxwellBlochField{N}
    E::ElectricField{N}
    P::PolarizationField{N}
    D::InversionField{N}
    function MaxwellBlochField(
                E::ElectricField{N},
                P::PolarizationField{N},
                D::InversionField{N}) where N
        return new{N}(E,P,D)
    end
end
function MaxwellBlochField(sim::Simulation)
    E = ElectricField(sim.x,zeros(ComplexF64,3length(sim),1))
    P = PolarizationField(sim.x,zeros(ComplexF64,3length(sim),1))
    D = InversionField(sim.x,zeros(Float64,length(sim)))
    G = AuxilliaryField(sim.x,zeros(ComplexF64,3length(sim),1))
    MaxwellBlochField(E,P,D,G)
end

function Base.show(io::IO,mbf::MaxwellBlochField)
    n = length(mbf.D)
    printstyled(io,"MaxwellBlochField",color=PRINTED_COLOR_DARK)
    print(io," (")
    printstyled(io,n,color=PRINTED_COLOR_NUMBER)
    print(io," points)")
end

Base.length(mb::MaxwellBlochField) = length(mb.E.pos)

struct MaxwellBloch{N,TM,TS}
    maxwell::TM
    source::TS
    fields::MaxwellBlochField{N}
    differentials::MaxwellBlochField{N}
end

function MaxwellBloch(sim::Simulation{1},ω::Real)
    sc=0
    m = maxwell(sim)
    m(ω)
    fields = MaxwellBlochField(sim)
    differentials = MaxwellBlochField(sim)
    return MaxwellBloch{1,typeof(m),typeof(sc)}(m,sc,fields,differentials)
end
# maxwell_bloch(args...;kwargs...) = MaxwellBoch(args..)

function Base.show(io::IO,mb::MaxwellBloch)
    n = length(mb.maxwell.sim)
    printstyled(io,"MaxwellBloch",color=PRINTED_COLOR_DARK)
    print(io," (")
    printstyled(io,n,color=PRINTED_COLOR_NUMBER)
    print(io," points)")
end

function (mb::MaxwellBloch)(dt::Number)
    # g=1
    # ħ=1
    # γparallel=1
    sim=mb.maxwell.sim
    # du[:] .= 0
    # u_to_fields!(mb,u)

    fields = mb.fields
    diffs = mb.differentials
    E, P, D, G = fields.E, fields.P, fields.D, fields.G
    dE, dP, dD, dG = diffs.E, diffs.P, diffs.D, diffs.G

    for i ∈ eachindex(dE.val)
        # k = mod1(i,length(mb.maxwell.sim))
        # d = mb.maxwell.sim.domain_indices[k]

        # if !isempty(sim.domains[1].χ)
            # γp = sim.domains[1].χ[1].γp
            # ωa = sim.domains[1].χ[1].ωa
        # else
            # γp = 0
            # ωa = 0
        # end

        dE[i] = -(mb.maxwell.αε[i,i])\G[i]# + 4π*(complex(γp,ωa))*P[i] - 4π*g^2*E[i]*D[i]/ħ)

        # dP.val[i] = 0#-complex(γp,ωa)*P[i] - g^2*E[i]*D[i]/ħ

    end

    mul!(dG.val,mb.maxwell.D²,E.val)
    # for i ∈ eachindex(dG.val) dG.val[i] /= mb.maxwell.αε[i,i] end

    for i ∈ eachindex(dE.val) E[i] += dt*dE[i] end
    for i ∈ eachindex(dG.val) G[i] += dt*dG[i] end

    # Ex, Ey, Ez = E.x, E.y, E.z
    # Px, Py, Pz = P.x, P.y, P.z
    # for i ∈ eachindex(dD.val)
        # if !isempty(sim.domains[1].χ)
            # D₀ = sim.domains[1].χ[1].D₀
        # else
            # D₀ = 0.0
        # end
        # dD.val[i] = γparallel*(D₀*sim.F[i]-D[i]) - 4*imag(Ex[i]*conj(Px[i])+Ey[i]*conj(Py[i])+Ez[i]*conj(Pz[i]))/ħ
    # end
    # fields_to_du!(du,mb)
    return nothing
end

end # module

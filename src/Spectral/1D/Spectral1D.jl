using Flatten
using SparseArrays

"""
    MaxwellNEP(::Simulation{1}; ky=0, kz=0) -> nep
"""
function MaxwellNEP(
            sim::Simulation{1,C,T},
            args...;
            ky::Number = 0,
            kz::Number = 0,
            kwargs...
            ) where {C,T}

    M = maxwell(sim; ky=ky, kz=kz)
    ka = 0
    f = map(f->(ω->f(ω,ka,ky,kz)),M.sim.self_energy.f)
    Σ1,Σ2,Σ3 = M.Σs
    αF, fχ = _compute_αχ(sim,args...)
    As = vcat([Σ1,Σ2,Σ3,sim.curlcurl(ky,kz),M.αε],αF)
    fs = [f[1],f[2],one,one,ω->-ω^2,map(ϝ->(ω->-ω^2*ϝ(ω)),fχ)...]
    nep = SPMF_NEP(As,fs;kwargs...)
    return MaxwellNEP(M,nep)
end


################################################################################
function _compute_αχ(sim::Simulation{1},args...)
    χs = flatten(map(d->d.χ,sim.domains),Common.Dispersions.AbstractDispersion)
    domain_inds = flatten(map(d->ntuple(i->d,length(sim.domains[d].χ)),ntuple(identity,length(sim.domains))))
    fχs = map(x->(ω->susceptability(x,ω)),χs)

    n = length(sim.x)
    Fs = map(d->sim.domains[d].pump,domain_inds)
    αFs = Array{SparseMatrixCSC{ComplexF64,Int},1}(undef,length(domain_inds)) # number of dispersions × number of domains
    for d ∈ eachindex(αFs)
        αFs[d] = spdiagm(0=>zeros(ComplexF64,3n))
        typeof(χs[d])<:TwoLevelSystem ? Common.Dispersions.TwoLevelSystems.H!(χs[d],args...) : nothing
        for i ∈ 1:3n
            j = mod1(i,n)
            if domain_inds[d]==sim.domain_indices[j]
                this_domain = sim.domain_indices[j]
                αFs[d].nzval[i] = Fs[this_domain](sim.x[j])*sim.α[1][j]*(isempty(args) ? 1 : χs[this_domain].hp1⁻¹[j])
            end
        end
    end
    return αFs,fχs
end








# @recipe function f(sim::Simulation,ψ::Array,inds=1:size(ψ,2); by=nothing, structure=false)
# 	@assert issubset(inds,1:size(ψ,2)) "indices $inds inconsistent with size(ψ,2)=$(size(ψ,2))"
# 	legend --> false
# 	aspect_ratio --> 1
# 	n = length(inds)
#     if isnothing(by)
#         if structure
#             layout --> (1+n,3)
# 			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/sqrt(9+n^2+1)
#         else
#             layout --> (n,3)
# 			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/sqrt(9+n^2)
#         end
#     else
#         if structure
#             layout --> (3+n)
# 			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/sqrt(9+n^2)
#         else
#             layout --> (n)
# 			markersize --> MARKERSIZE_SCALE/sqrt(size(ψ,1))/n
#         end
#     end
# 	if structure
# 		if isnothing(by)
# 			@series (sim,by)
# 		else
# 			@series (sim,real)
# 		end
# 	end
# 	for i ∈ 1:n
# 		if isnothing(by)
# 			bys = [:real, :imag, :abs2]
# 			for j ∈ eachindex(bys)
# 				if structure
# 					@series begin
# 						subplot --> 3i+j
# 						(sim,ψ,inds[i],bys[j])
# 					end
# 				else
# 					@series begin
# 						subplot --> 3(i-1)+j
# 						(sim,ψ,inds[i],bys[j])
# 					end
# 				end
# 			end
# 		else
# 			@series begin
# 				subplot --> i
# 				colorbar --> false
# 				(sim,ψ,inds[i],by)
# 			end
# 		end
# 	end
# end
# @recipe function f(sim::Simulation,ψ::Array,ind::Int,by::Union{Symbol,Function})
# 	markershape --> :rect
# 	markerstrokealpha --> 0
# 	seriestype --> :scatter
# 	q = quantile(abs.(ψ[:,ind]),PLOT_QUANTILE)
# 	if by ∈ [:real, :Real, :re, :Re, real]
# 		@series begin
# 			title --> "Re psi"
# 			markercolor --> :diverging
# 			colorbar --> false
# 			clim --> (-q,q).*PLOT_SCALE_FUDGE
# 			z = real(ψ[:,ind])
# 			marker_z --> z
# 			(sim.x[sim.interior],sim.y[sim.interior])
# 		end
# 	elseif by ∈ [:imag, :Imag, :im, :Im, imag]
# 		@series begin
# 			title --> "Im psi"
# 			markercolor --> :diverging
# 			colorbar --> false
# 			clim --> (-q,q).*PLOT_SCALE_FUDGE
# 			z = imag(ψ[:,ind])
# 			marker_z --> z
# 			(sim.x[sim.interior],sim.y[sim.interior])
# 		end
# 	elseif by ∈ [:abs2, :Abs2, abs2]
# 		@series begin
# 			title --> "|psi|²"
# 			markercolor --> :sequential
# 			colorbar --> false
# 			clim --> (0,q^2).*PLOT_SCALE_FUDGE
# 			z = abs2.(ψ[:,ind])
# 			marker_z --> z
# 			(sim.x[sim.interior],sim.y[sim.interior])
# 		end
# 	elseif by ∈ [:abs, :Abs, abs]
# 		@series begin
#             title --> "|psi|"
# 			markercolor --> :sequential
# 			colorbar --> false
# 			z = abs.(ψ[:,ind])
# 			clim --> (0,q).*PLOT_SCALE_FUDGE
# 			marker_z --> z
# 			(sim.x[sim.interior],sim.y[sim.interior])
# 		end
# 	end
# end

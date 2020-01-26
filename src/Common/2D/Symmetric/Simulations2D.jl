# TODO: fix and test smoothing
module Simulations2DSymmetric


using ..Curlcurls
using ..Dispersions
using ..Domains
using ..Laplacians
using ..Lattices
using ..Points
using ..SelfEnergies
using ..Shapes
using ..VectorFields
# using LinearAlgebra
using SparseArrays
using Statistics

import ..AbstractDomain
import ..Symmetric, ..Unsymmetric
import ..Simulation
import ..smooth_dielectric!
import ..smooth_pump!
import ..which_domains
import ..simulation_dielectric!
import ..simulation_pump!

function Simulation(
			ω₀::Real,
			lattice_domains::Tuple{LatticeDomain{2,Symmetric,Cartesian}},
			nondispersive_domains::NTuple{L,NondispersiveDomain{2}},
			dispersive_domains::NTuple{M,DispersiveDomain{2}};
			k₃₀::Real=0
			) where {L,M}

	lattice_domain = lattice_domains[1] # there is only one lattice by construction, and it's Cartesian
	smoothed = Ref(false)

	lattice = lattice_domain.lattice
	indices = lattice_domain.indices
	imin, imax = extrema(map(ld->ld[1],indices))
	jmin, jmax = extrema(map(ld->ld[2],indices))
	x = [lattice[i,0] for i ∈ imin:imax]
	y = [lattice[0,j] for j ∈ jmin:jmax]

	lattice_domain_indices = fill(1,length(lattice_domain.x))
	nondispersive_domain_indices = which_domains(nondispersive_domains,lattice_domain.x)
	dispersive_domain_indices = which_domains(dispersive_domains,lattice_domain.x)

	# populate dielectric
	ε = Vector{ComplexF64}(undef,length(lattice_domain.x))
	if isempty(nondispersive_domains)
		for i ∈ eachindex(ε) ε[i] = lattice_domain.ε end
	else
		dielectrics = map(n -> getfield(n,:dielectric), nondispersive_domains)
		for i ∈ eachindex(ε)
			d = nondispersive_domain_indices[i]
			ε[i] = iszero(d) ? lattice_domain.ε : dielectrics[d](lattice_domain.x[i])
		end
	end

	# populate pump and dispersive susceptability
	F = Vector{Float64}(undef,length(lattice_domain.x))
	χ = Vector{AbstractDispersion}(undef,length(lattice_domain.x))
	Fs = Vector{Vector{Float64}}(undef,length(dispersive_domains))
	for i ∈ eachindex(Fs) Fs[i] = zeros(Float64,length(x)) end
	if isempty(dispersive_domains)
		for i ∈ eachindex(F) F[i] = 0 end
		for i ∈ eachindex(F) χ[i] = NoDispersion() end
	else
		pumps = map(n -> n.pump, dispersive_domains)
		χs = map(n -> n.χ, dispersive_domains)
		for i ∈ eachindex(F)
			d = dispersive_domain_indices[i]
			if iszero(d)
				F[i] = 0
				χ[i] = NoDispersion()
			else
				F[i] = pumps[d](lattice_domain.x[i])
				χ[i] = χs[d]
				Fs[d][i] = pumps[d](x[i])
			end
		end
	end

	# generate boundary layers
	σx, _ = boundary_layer(lattice_domain, x)
	_, σy = boundary_layer(lattice_domain, y)
	αx, αy = 1 .+ 1im*σx/sqrt(ω₀^2-k₃₀^2), 1 .+ 1im*σy/sqrt(ω₀^2-k₃₀^2)

	# generate half sites and boundary layers
	x_half, y_half = _generate_half_xy(lattice_domain, x)
	σx_half, _ = boundary_layer(lattice_domain, x_half)
	_, σy_half = boundary_layer(lattice_domain, y_half)
	αx_half, αy_half = 1 .+ 1im*σx_half/sqrt(ω₀^2-k₃₀^2), 1 .+ 1im*σy_half/sqrt(ω₀^2-k₃₀^2)

	laplacian = Laplacian{Symmetric}(lattice_domain.lattice, αx, αy, αx_half, αy_half)
	# curlcurl = Curlcurl{Symmetric}(lattice_domain.lattice,α[1],α_half[1],nnm,nnp,indices,interior,surface)

	Σ = SelfEnergy{Symmetric}(lattice_domain, αx_half, αy_half)

	return Simulation{2,Symmetric,ComplexF64,typeof(lattice_domains),typeof(nondispersive_domains),typeof(dispersive_domains),typeof(Σ)}(
		lattice_domains,
		nondispersive_domains,
		dispersive_domains,
		lattice_domain.x,
		lattice_domain_indices,
		nondispersive_domain_indices,
		dispersive_domain_indices,
		ε,
		F,
		χ,
		Fs,
		laplacian,
		# curlcurl,
		Σ,
		(αx, αy),
		(σx, σy),
		(x_half, y_half),
		(σx_half, σy_half),
		ω₀,
		NaN,
		NaN,
		k₃₀,
		smoothed)
end

################################################################################
# SIMULATION{2} building utilities

# used in building Simulation{2}
function _generate_half_xy(lattice_domain::LatticeDomain{2,Symmetric,Cartesian}, x::Vector{TP}) where TP<:Point{2}
	lattice = lattice_domain.lattice
	indices = lattice_domain.indices
	imin, imax = extrema(map(ld->ld[1],indices))
	jmin, jmax = extrema(map(ld->ld[2],indices))
	half_x = [lattice[i,0] for i ∈ ((imin:(imax+1)).-1/2)]
	half_y = [lattice[0,j] for j ∈ ((jmin:(jmax+1)).-1/2)]
	return half_x, half_y
end

@inline function boundary_layer(domain::LatticeDomain{2,Symmetric,Cartesian}, x::Vector{Point{2,Cartesian}})
	σx = domain.boundary.bls[1].(x) .+ domain.boundary.bls[2].(x)
	σy = domain.boundary.bls[3].(x) .+ domain.boundary.bls[4].(x)
	return σx, σy
end


################################################################################
# 2-d smoothing functions

function smooth_dielectric!(sim::Simulation{2,Symmetric}, num_sub_pixel::Integer=NUM_SUBPIXELS)
    indices = sim.nondispersive_domain_indices
    lattice = sim.lattice
    X = Matrix{Point{2,Cartesian}}(undef,num_sub_pixel,num_sub_pixel)
    E = Vector{ComplexF64}(undef,num_sub_pixel)
	r = LinRange(-.5,.5,num_sub_pixel)
    for i ∈ eachindex(indices)
		if 1<i<length(indices)
			idx, = Tuple(sim.lattice_domain.indices[i])
			if !(indices[i] == indices[i-1] == indices[i+1])
				for j ∈ eachindex(r) X[j] = lattice[idx+r[j]] end
		        simulation_dielectric!(E,sim,X)
				sim.ε[i] = mean(E)
			end
		else
			# placeholder for periodic lattice domains
		end
    end
    return nothing
end

function smooth_pump!(sim::Simulation{2}, num_sub_pixel::Integer = NUM_SUBPIXELS)
	indices = sim.dispersive_domain_indices
    lattice = sim.lattice
    X = Vector{Point{1}}(undef,num_sub_pixel)
    F = Vector{ComplexF64}(undef,num_sub_pixel)
	r = LinRange(-.5,.5,num_sub_pixel)
    for i ∈ eachindex(indices)
		if 1<i<length(indices)
			idx, = Tuple(sim.lattice_domain.indices[i])
			if !(indices[i] == indices[i-1] == indices[i+1])
				for j ∈ eachindex(r) X[j] = lattice[idx+r[j]] end
		        simulation_pump!(F,sim,X)
				sim.F[i] = mean(F)
			end
		end
    end
    return nothing
end

################################################################################
# extras

function Base.getproperty(sim::Simulation{2}, sym::Symbol)
	if sym==:lattice_domain
		return getfield(sim,:lattice_domains)[1]
	elseif Base.sym_in(sym,(:k10,:k1₀))
		return getfield(sim,:k₁₀)
	elseif Base.sym_in(sym,(:k20,:k2₀))
		return getfield(sim,:k₂₀)
	elseif Base.sym_in(sym,(:k30,:k3₀))
		return getfield(sim,:k₃₀)
	elseif Base.sym_in(sym,(:shape,:boundary,:lattice))
		return getproperty(getfield(sim,:lattice_domains)[1],sym)
	elseif Base.sym_in(sym,propertynames(getproperty(getfield(sim,:lattice_domains)[1],:lattice)))
		return getproperty(getfield(sim,:lattice_domains)[1],sym)
	elseif Base.sym_in(sym,(:lat,:lattice,:Lat,:Lattice))
		return getfield(sim,:lattice_domains)[1].lattice
	else
		return getfield(sim,sym)
	end
end

function Base.propertynames(sim::Simulation{2,Symmetric}, private=false)
	if private
		return fieldnames(Simulation)
	else
		return (:shape, :boundary, :lattice, :lattice_domain, :lattice_domains, :nondispersive_domains, :dispersive_domains, :ε, :F, :x, :χ, :ω₀, :k1₀, :k2₀, :k3₀, :lattice, propertynames(sim.lattice)...)
	end
end

################################################################################
# Pretty Printing

import ...PRINTED_COLOR_GOOD
import ...PRINTED_COLOR_WARN
import ...PRINTED_COLOR_DARK

function Base.show(io::IO,sim::Simulation{2,CLASS}) where CLASS
	if sim.smoothed[]
		printstyled(io,"Smoothed ",color=PRINTED_COLOR_GOOD)
	else
		printstyled(io,"Unsmoothed ",color=PRINTED_COLOR_WARN)
	end
	print(io,"2D ")
	CLASS<:Symmetric ? print(io,"Symmetric") : nothing
	CLASS<:Unsymmetric ? print(io,"Unsymmetric") : nothing
	printstyled(io," Simulation\n",color=PRINTED_COLOR_DARK)
    print(io,"\tn sites: ")
	printstyled(io,length(sim),"\n",color=:light_cyan)
	println(io,"\t====================")
	println(IOContext(io,:tabbed2=>true),sim.boundary)
	println(io)
	println(IOContext(io,:tabbed2=>true),sim.lattice)
	println(io)
	domains = (sim.nondispersive_domains...,sim.dispersive_domains...)
	for d ∈ eachindex(domains)
		print(io,"\t\tDomain ", d, " ")
		print(IOContext(io,:tabbed2=>true),domains[d])
		d < length(domains) ? println(io) : nothing
	end
end

end

using .Simulations2DSymmetric

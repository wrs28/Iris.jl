# TODO: add waveguide sites for matched conditions
"""
	module Simulations

Package data used to define simulation and generate the relevant differential operators.
Exports are `Simulation`, `smooth!`, `update_dielectric!`, `update_pump`
"""
module Simulations

export Simulation
export smooth!
export update_dielectric!
export update_pump!
export Unsymmetric
export Symmetric
export Hermitian

files = (
	"1D/Simulations1D.jl",
	# "2D/Simulations2D.jl",
	# "3D/Simulations3D.jl"
	)

using ..Curlcurls
using ..Dispersions
using ..Domains
using ..Laplacians
using ..Lattices
using ..Points
using ..SelfEnergies
using ..Shapes
using ..VectorFields
using LinearAlgebra
using SparseArrays
using Statistics

import ..AbstractDomain

"""
	Simulation(ω₀, domains...; [k₂₀, k₃₀, k₁₀]) -> simulation

Build a simulation object which can be passed to the various solvers to find
resonance frequencies, etc.

There are three types of `domains`: `LatticeDomain`, `NondispersiveDomain`,
`DispersiveDomain`, which can be provided in any order, though at least one
`LatticeDomain` must be provided.
Where domains of the same type overlap, the one that appears first takes priority.
"""
struct Simulation{N,CLASS,T,TLDOM,TNDOM,TDDOM,TSE}
	lattice_domains::TLDOM
	nondispersive_domains::TNDOM
	dispersive_domains::TDDOM

	x::Vector{Point{N}}
	lattice_domain_indices::Vector{Int}
	nondispersive_domain_indices::Vector{Int}
	dispersive_domain_indices::Vector{Int}
	ε::Vector{ComplexF64}
	F::Vector{Float64}
	χ::Vector{AbstractDispersion}

	laplacian::Laplacian{N}
	curlcurl::Curlcurl{N}
	self_energy::TSE
	α::NTuple{N,Vector{ComplexF64}}
	σ::NTuple{N,Vector{ComplexF64}}
	x_half::NTuple{N,Vector{Point{N}}}
	σ_half::NTuple{N,Vector{ComplexF64}}
	ω₀::Float64
	k₁₀::Float64
	k₂₀::Float64
	k₃₀::Float64
	smoothed::Base.RefValue{Bool}
end

foreach(include,files)

Base.length(sim::Simulation) = length(sim.x)
Base.ndims(sim::Simulation{N}) where N = N

# Base.conj(sim::Simulation) = Simulation(sim.ω₀,map(conj,sim.domains)...)

import ..NUM_SUBPIXELS

"""
	smooth!(sim,[num_sub_pixel])

Use sub-pixel sampling along domain walls to smooth the dieletric and pump profiles.

`num_sub_pixel` defaults to `NUM_SUB_PIXEL` in module `Defaults`, which is currently set to 5.
"""
function smooth!(sim::Simulation, num_sub_pixel::Integer = NUM_SUBPIXELS)
    smooth_dielectric!(sim, num_sub_pixel)
	smooth_pump!(sim, num_sub_pixel)
	sim.smoothed[] = true
    return nothing
end

"""
	update!(sim)

update simulation `sim` by recomputing the dielectric and pump functions
"""
function VectorFields.update!(sim::Simulation)
	update_dielectric!(sim)
	update_pump!(sim)
	return nothing
end

"""
	update_dielectric!(sim,[ε])

update simulation `sim` by recomputing the dielectric function or by giving it an explicit array `ε`
"""
function update_dielectric!(sim::Simulation)
	simulation_dielectric!(sim.ε, sim, sim.x)
	sim.smoothed[] = false
	return nothing
end
function update_dielectric!(sim::Simulation,ε::AbstractVecOrMat)
	copyto!(sim.ε,ε)
	sim.smoothed[] = false
	return nothing
end

"""
	update_pump!(sim,[F])

update simulation `sim` by recomputing the pump function or by giving it an explicit array `F`
"""
function update_pump!(sim::Simulation)
	simulation_pump!(sim.F, sim, sim.x)
	sim.smoothed[] = false
	return nothing
end
function update_pump!(sim::Simulation,F::AbstractVecOrMat)
	copyto!(sim.F, F)
	sim.smoothed[] = false
	return nothing
end


# the dielectric as a function of position (x,y)
# function simulation_dielectric(sim::Simulation, x::Vector)
	# ε = Vector{ComplexF64}(undef,length(x))
	# simulation_dielectric!(ε, sim , x)
	# return ε
# end
function simulation_dielectric!(ε::Vector, sim::Simulation, x::Vector, backgrounds=map(d->d.ε, sim.lattice_domains), dielectrics=map(d->d.dielectric, sim.nondispersive_domains))
    for i ∈ eachindex(x) ε[i] = simulation_dielectric(sim, x[i], backgrounds, dielectrics) end
    return nothing
end
@inline function simulation_dielectric(sim::Simulation, x::Point, backgrounds=map(d->d.ε, sim.lattice_domains), dielectrics=map(d->d.dielectric, sim.nondispersive_domains))
    lattice_index = which_domain(sim.lattice_domains, x)
	nondispersive_index = which_domain(sim.nondispersive_domains, x)
	if nondispersive_index==0
		return backgrounds[lattice_index]
	else
	    return dielectrics[nondispersive_index](x)
	end
end

# the dielectric as a function of position (x,y)
# function simulation_pump(sim::Simulation, x::Vector)
# 	F = Vector{Float64}(undef,length(x))
# 	simulation_pump!(F, sim , x)
# 	return F
# end
function simulation_pump!(F::Vector, sim::Simulation, x::Vector, pumps=map(d->d.pump, sim.dispersive_domains))
    for i ∈ eachindex(x) F[i] = simulation_pump(sim, x[i], pumps) end
    return nothing
end
@inline function simulation_pump(sim::Simulation, x::Point, pumps=map(d->d.pump, sim.dispersive_domains))
	dispersive_index = which_domain(sim.dispersive_domains, x)
	if dispersive_index==0
		return 0.0
	else
	    return pumps[dispersive_index](x)
	end
end

################################################################################
# utilities

# identify which domain is associated with a given spatial position (x)
function which_domains(domains, x::Vector)
    indices = Vector{Int}(undef, length(x))
	shapes = map(z->z.shape, domains)
    which_domains!(indices, domains, shapes, x)
    return indices
end
function which_domains!(indices::Vector{Int}, domains, shapes, x)
    @inbounds for i ∈ eachindex(x) indices[i] = which_domain(domains, shapes, x[i]) end
    return nothing
end
# above is for arrays of points, below is for single points
function which_domain(domains,x::Point)
	shapes = map(z->z.shape,domains)
    return which_domain(domains, shapes, x)
end
@inline function which_domain(domains, shapes, x::Point)
	index = findfirst(z->z(x),shapes)
    isnothing(index) ? index=0 : nothing
    iszero(index) ? nothing : (isvoid(domains[index]) ? index=0 : nothing)
    return index
end

isvoid(domain::AbstractDomain) = domain.type==:Void

function boundary_layer!(Σ, domains, x::Vector{Point{N}}, domain_indices) where N
	for j ∈ eachindex(Σ) for i ∈ eachindex(x)
		Σ[j][i] = boundary_layer(domains[domain_indices[i]],x[i],j)
	end end
	return nothing
end












# take individually defined domains and concatenate all relevant data
function append_domains!(domain_index,nnm,nnp,name,x,indices,ε,F,interior,bulk,surface,corner,domains)
	count = 1
	for d ∈ domains
		append!(domain_index,fill(count,length(d.x)))
		for i ∈ eachindex(d.nnm)
			append!(nnm[i],d.nnm[i].+length(nnm[i]))
			append!(nnp[i],d.nnp[i].+length(nnp[i]))
		end
		append!(name,fill(d.name,length(d.x)))
		append!(x,d.x)
		append!(indices,d.indices)
		append!(ε,d.ε)
		append!(F,d.F)
		append!(interior,d.interior)
		append!(bulk,d.bulk)
		append!(surface,d.surface)
		append!(corner,d.corner)
		count += 1
	end
	return nothing
end

# identify which sites are to be removed by determining which sites belong to the topmost domain.
get_sites_to_remove(domains,dom,x) = which_domains(domains,x).!==dom

# remove the pre-computed nearest neighbors, which makes the bulk nn calculation very fast
function trim_nn!(nnm,nnp,removed,bulk)
	cr = cumsum(removed)
	for j ∈ eachindex(nnm)
		for i ∈ eachindex(nnm[j])
			if bulk[i]
				nnm[j][i] -= cr[nnm[j][i]]
				nnp[j][i] -= cr[nnp[j][i]]
			end
		end
	end
end

# remove the sites determined by get_sites_to_remove
remove_sites!(sites::Tuple,removed) = map(z->remove_sites!(z,removed),sites)
remove_sites!(site::NTuple{N,Array},removed) where N = map(z->deleteat!(z,findall(removed)),site)
remove_sites!(site,removed) = deleteat!(site,findall(removed))


################################################################################
#Pretty Printing

import ..PRINTED_COLOR_GOOD
import ..PRINTED_COLOR_WARN
import ..PRINTED_COLOR_DARK

function Base.show(io::IO,sim::Simulation{N,CLASS}) where {N,CLASS}
	suffix = length(sim.domains)>1 ? "s" : ""
	if sim.smoothed[]
		printstyled(io,"Smoothed ",color=PRINTED_COLOR_GOOD)
		print(io,"$(N)D ", CLASS)
		printstyled(io," Simulation",color=PRINTED_COLOR_DARK)
		println(io, " with ", length(sim.domains)," domain",suffix)
	else
		printstyled(io,"Unsmoothed ",color=PRINTED_COLOR_WARN)
		print(io,N, "D ", CLASS)
		printstyled(io," Simulation",color=PRINTED_COLOR_DARK)
		println(io, " with ", length(sim.domains)," domain",suffix)
	end
    print(io,"\tn sites: ")
	printstyled(io,length(sim),"\n",color=:light_cyan)
	println(io,"\t\t====================")
	for d ∈ eachindex(sim.domains)
		printstyled(io,"\t\tDomain ",color=PRINTED_COLOR_DARK)
		print(io, d," (",sim.domains[d].type,"): ")
		print(io,sim.domains[d].name)
		d < length(sim.domains) ? println(io) : nothing
	end
end


end # module

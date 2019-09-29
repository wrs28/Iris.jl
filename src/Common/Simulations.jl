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
export update_sim
export Unsymmetric
export Symmetric
export Hermitian

using ...Defaults
using ..Curlcurls
using ..Domains
using ..Lattices
using ..Points
using ..SelfEnergies
using ..Shapes
using LinearAlgebra
# using RecipesBase
using SparseArrays

struct Unsymmetric end

"""
	struct Simulation

	Simulation(domains...; verbose=true) -> sim

	(::Simulation)(;kwargs...) -> new_sim

Takes all the given `domains` and builds a simulation object which can be passed
to the various solvers to find resonance frequencies, etc. The where domains overlap,
the one that appears first in the arguments of `Simulation` takes priority.

Contains all the information that defines a simulation (except for frequency).
It contains several kinds of fields:
 * vectors of data about each site, such as `x`, `y`, `ε`, `F`
 * boolean vectors that give geometric properites of each site, such as whether it is an `interior` point
 * vectors of weights and indices giving the coupling between sites
 * vectors of weights and indices giving the self-couplings
 * the `tessellation` used to guide interpolation

 Constructing a `Simulation` is sped up if multiple threads are running,
 i.e. by declaring the environment variable "JULIA\\_NUM\\_THREADS" *before*
 starting Julia.
"""
struct Simulation{N,C,T,TDOM,M,TF}
	smoothed::Base.RefValue{Bool}
	ε::Array{ComplexF64,1}
	F::Array{Float64,1}
	curlcurl::Curlcurl{N}
	self_energy::SelfEnergy{N,M,TF}
	domains::TDOM
	domain_indices::Array{Int,1}
	x::Array{Point{N},1}
	x_half::NTuple{N,Array{Point{N},1}}
	α::NTuple{N,Array{ComplexF64,1}}
	σ::NTuple{N,Array{ComplexF64,1}}
	σ_half::NTuple{N,Array{ComplexF64,1}}
	name::Array{Symbol,1}
	ω₀::Float64
	k₁₀::Float64
	k₂₀::Float64
	k₃₀::Float64
end

include("1D/Simulations.jl")
# include("2D/Simulations.jl")
# include("3D/Simulations.jl")

function Base.show(io::IO,sim::Simulation{N,CLASS}) where {N,CLASS}
	suffix = length(sim.domains)>1 ? "s" : ""
	if sim.smoothed[]
        println(io,"Smoothed ", N, "D ", CLASS, " Simulation with ", length(sim.domains) ," domain", suffix)
	else
		println(io,"Unsmoothed ", N, "D ", CLASS, " Simulation with ", length(sim.domains)," domain",suffix)
	end
    println(io,"\tsites: ",length(sim.x))
	println(io,"\t\t====================")
	for d ∈ eachindex(sim.domains)
		println(io,"\t\tDomain ", d,": ",sim.domains[d].type)
	end
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

# identify which domain is associated with a given spatial position (x)
function which_domains(domains,x)
    doms = Array{Int}(undef,length(x))
	shapes = map(z->z.boundary.shape,domains)
    which_domains!(doms,domains,shapes,x)
    return doms
end
function which_domains!(doms::Array{Int},domains,shapes,x)
    for i ∈ eachindex(x)
		doms[i] = which_domain(domains,shapes,x[i])
    end
    return nothing
end
function which_domain(domains,x)
	bnd = map(z->z.boundary,domains)
	shapes = map(z->z.boundary.shape,domains)
    return which_domain(domains,shapes,x)
end
@inline function which_domain(domains,shapes,x::Point)
	dom = findfirst(map(z->z(x),shapes))
    isnothing(dom) ? dom=0 : nothing
    iszero(dom) ? nothing : (isvoid(domains[dom]) ? dom=0 : nothing)
    return dom
end

isvoid(domain::Domain) = domain.type==:Void

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

"""
	update_dielectric!(sim,[ε])

update simulation `sim` by recomputing the dielectric function or by giving it an explicit array `ε`
"""
function update_dielectric!(sim::Simulation,ε::Array)
	sim.ε[:] = ε
	sim.smoothed[] = false
	return nothing
end
update_dielectric!(sim::Simulation) = update_dielectric!(sim,simulation_dielectric.(Ref(sim),sim.x))
"""
	update_pump!(sim,[F])

update simulation `sim` by recomputing the pump function or by giving it an explicit array `F`
"""
function update_pump!(sim::Simulation,F::Array)
	sim.F[:] = F
	sim.smoothed[] = false
	return nothing
end
update_pump!(sim::Simulation) = update_pump!(sim,simulation_pump.(Ref(sim),sim.x))


"""
	update_sim!(sim)

update simulation `sim` by recomputing the dielectric and pump functions
"""
function update_sim!(sim::Simulation)
	update_dielectric!(sim)
	update_pump!(sim)
	return nothing
end


"""
	smooth!(sim,[num_sub_pixel])

Use sub-pixel sampling along domain walls to smooth the dieletric and pump profiles.

`num_sub_pixel` defaults to `NUM_SUB_PIXEL` in module `Defaults`, which is currently set to 5.
"""
smooth!

# used in the smoothing routing smooth!, it computes the dielectric as a function of position (x,y)
function simulation_dielectric!(ε,domains,x)
    for i ∈ eachindex(x)
        ε[i] = simulation_dielectric(domains,x[i])
    end
    return nothing
end
function simulation_dielectric(domains,x)
    dom = which_domain(domains,x)
    ε = dom==0 ? complex(NaN,NaN) : domains[dom].dielectric(x)
    return ε
end
simulation_dielectric(sim::Simulation,x) = simulation_dielectric(sim.domains,x)

# used in the smoothing routing smooth!, it computes the pump as a function of position (x,y)
function simulation_pump!(F,domains,x)
    for i ∈ eachindex(x)
        F[i] = simulation_pump(domains,x[i])
    end
    return nothing
end
function simulation_pump(domains,x)
    dom = which_domain(domains,x)
    F = dom==0 ? 0.0 : domains[dom].pump(x)
    return F
end
simulation_pump(sim::Simulation,x) = simulation_pump(sim.domains,x)



function boundary_layer(domains,x::Array{Point{N}},domain_index) where N
	Σ = ntuple(i->Array{ComplexF64}(undef,length(x)),N)
	boundary_layer!(Σ,domains,x,domain_index)
	return Σ
end
function boundary_layer!(Σ,domains,x::Array{Point{N}},domain_index) where N
	for j ∈ eachindex(Σ)
		for i ∈ eachindex(x)
			Σ[j][i] = boundary_layer(domains[domain_index[i]],x[i],j)
		end
	end
	return nothing
end

end # module

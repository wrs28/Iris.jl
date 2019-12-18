# TODO: add waveguide sites for matched conditions

"""
	module Simulations

Package data used to define simulation and generate the relevant differential operators.
Exports are `Simulation`, `smooth!`, `update_dielectric!`, `update_pump`
"""
module Simulations

const N_NEAREST_SITES_FOR_TESSELLATION = 15

export Simulation,
smooth!,
update_dielectric!,
update_pump!,
update_sim!

using ...Defaults
using ..Boundaries
using ..Domains
using ..Lattices
using ..Shapes
using ..Tessellations
using Combinatorics
using LinearAlgebra
using ProgressMeter
using RecipesBase
using SparseArrays
using Statistics
using VoronoiDelaunay


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
struct Simulation
	domains::Array{Domain,1}
	nsites::Int
	n::Int

	# site boolean properties, one for each site
	interior::BitArray{1}
	exterior::BitArray{1}
	corner::BitArray{1}
	bulk::BitArray{1}
	surface::BitArray{1}
	lattice_wall::BitArray{1}
	domain_wall::BitArray{1}

	# site properties, indices
	domain_index::Array{Int,1}
	nnxm::Array{Int,1}
	nnxp::Array{Int,1}
	nnym::Array{Int,1}
	nnyp::Array{Int,1}
	ij::Array{CartesianIndex{2},1}
	exterior_from::Array{Int,1}
	exterior_side::Array{Int,1}
	not_this_domains_exterior::Array{Array{Int,1},1}
	exterior_inds::Array{Int,1}


	# site properties, other
	x::Array{Float64,1}
	y::Array{Float64,1}
	r::Array{Float64,1}
	ε::Array{ComplexF64,1}
    F::Array{Float64,1}
	name::Array{Symbol,1}
	distance::Array{Array{Float64,1},1}
	normal::Array{Array{Array{Float64,1},1},1}
	tangent::Array{Array{Array{Float64,1},1},1}
	side::Array{Array{Int,1},1}

	# tessellation::Tessellation
	smoothed::Base.RefValue{Bool}

	# links between sites
	link_start::Array{Int,1}
	link_stop::Array{Int,1}
	link_weight::Array{Float64,1}
	link_half_x::Array{Float64,1}
	link_half_y::Array{Float64,1}
	link_half_r::Array{Float64,1}
	link_x_bool::BitArray{1}
	link_y_bool::BitArray{1}
	link_r_bool::BitArray{1}

	# self-energy terms
	self_index::Array{Int,1}
	self_weight::Array{Float64,1}
	self_half_x::Array{Float64,1}
	self_half_y::Array{Float64,1}
	self_half_r::Array{Float64,1}
	self_x_bool::BitArray{1}
	self_y_bool::BitArray{1}
	self_r_bool::BitArray{1}

	Simulation(sim::Simulation) = Simulation(sim.domains;kwargs...)
    Simulation(domains...; kwargs...) = Simulation(domains; kwargs...)
    function Simulation(domains::Tuple; verbose::Bool=false)

		domain_index = Int[]
		nnxm, nnxp, nnym, nnyp = Int[], Int[], Int[], Int[]
		name = Symbol[]
		x, y, r, θ = Float64[], Float64[], Float64[], Float64[]
		ij = CartesianIndex{2}[]
		ε, F = ComplexF64[], Float64[]
		interior, surface, bulk = BitArray([]), BitArray([]), BitArray([])
		corner = BitArray([])

		verbose ? println("\nBuilding Simulation...") : nothing
		verbose ? println("\tCombining $(length(domains)) domains") : nothing
			append_domains!(domain_index,nnxm,nnxp,nnym,nnyp,name,x,y,r,θ,ij,ε,F,interior,bulk,surface,corner,domains)

			removed = get_sites_to_remove(domains,domain_index,x,y)

			trim_nn!(nnxm,nnxp,nnym,nnyp,removed,bulk)

			remove_sites!((domain_index,name,x,y,r,θ,ij,ε,F,interior,bulk,surface,corner,nnxm,nnxp,nnym,nnyp),removed)

		verbose ? println("\tComputing $(sum(surface)) normals and extending domain") : nothing
			lattice_wall, domain_wall = generate_lattice_and_domain_wall!(interior,bulk,surface,corner,domains,domain_index,x,y,ij)

			exterior_x, exterior_y, exterior_from, normal, tangent, distance, side, exterior_side = generate_exterior!(corner,domains,domain_index,x,y,surface,verbose)

			exterior = falses(length(x))

			append_exterior!(domain_index,nnxm,nnxp,nnym,nnyp,name,x,y,r,θ,ij,ε,F,interior,bulk,surface,corner,
				exterior,domain_wall,lattice_wall,normal,tangent,distance,side,exterior_from,exterior_side,exterior_x,exterior_y)

			not_this_domains_exterior = Array{Array{Int,1}}(undef,length(domains))#falses(length(x),length(domains))
			for d ∈ eachindex(domains)
				not_this_domains_exterior[d] = sum(interior) .+ findall(!isequal(d),domain_index[exterior_from[exterior]])
			end
			exterior_inds = findall(exterior)

		# verbose ? println("\tTessellating $(length(x)) points") : nothing
			# tessellation = Tessellation(x,y)

		bulk_link_size = 4sum(bulk)
		lattice_wall_link_size = 4NNN_BULK*sum(lattice_wall)
		surface_link_size = (NNN_BULK+2NNN_SURFACE)*sum(surface)
		surface_self_size = sum(surface)
		corner_link_size = 2NNN_CORNER*sum(corner)
		corner_self_size = 2sum(corner)

		link_start = Array{Int}(undef,bulk_link_size+lattice_wall_link_size+surface_link_size+corner_link_size)
		link_stop = Array{Int}(undef,length(link_start))
		link_weight = Array{Float64}(undef,length(link_start))
		link_half_x = Array{Float64}(undef,length(link_start))
		link_half_y = Array{Float64}(undef,length(link_start))
		link_half_r = Array{Float64}(undef,length(link_start))
		link_x_bool = falses(length(link_start))
		link_y_bool = falses(length(link_start))
		link_r_bool = falses(length(link_start))

		self_index = Array{Int}(undef,surface_self_size+corner_self_size)
		self_weight = Array{Float64}(undef,length(self_index))
		self_half_x = Array{Float64}(undef,length(self_index))
		self_half_y = Array{Float64}(undef,length(self_index))
		self_half_r = Array{Float64}(undef,length(self_index))
		self_x_bool = falses(length(self_index))
		self_y_bool = falses(length(self_index))
		self_r_bool = falses(length(self_index))

		link_fiduciary = falses(length(link_start))
		self_fiduciary = falses(length(self_index))

		sim = new([domains...],length(x),sum(interior),interior,exterior,corner,bulk,
			surface,lattice_wall,domain_wall,
			domain_index,nnxm,nnxp,nnym,nnyp,ij,
			exterior_from,exterior_side,not_this_domains_exterior,exterior_inds,
			x,y,r,ε,F,name,distance,normal,tangent,side,
			# tessellation,
			Ref(false),
			link_start,link_stop,link_weight,
			link_half_x,link_half_y,link_half_r,link_x_bool,link_y_bool,link_r_bool,
			self_index,self_weight,
			self_half_x,self_half_y,self_half_r,self_x_bool,self_y_bool,self_r_bool)

		verbose ? println("\tLinking $(sum(sim.bulk)) bulk points") : nothing
			link_bulk!(link_fiduciary,sim,verbose)

		verbose ? println("\tLinking $(sum(sim.lattice_wall)) lattice mismatch points") : nothing
			link_lattice_wall!(link_fiduciary,sim,verbose)

		verbose ? println("\tLinking $(sum(sim.surface)) surface points") : nothing
			link_surface!(link_fiduciary,self_fiduciary,sim,verbose)

		verbose ? println("\tLinking $(sum(sim.corner)) corner points") : nothing
			link_corner!(link_fiduciary,self_fiduciary,sim,verbose)

		verbose ? println("\tCleaning up $(sum(link_fiduciary)+sum(self_fiduciary)) links") : nothing
			trim_link!(sim,link_fiduciary)

			trim_self!(sim,self_fiduciary)

		verbose ? println("\tApplying boundary conditions") : nothing
			apply_local_bc!(sim)

			apply_nonlocal_bc!(sim)

		verbose ? println("\tDone!") : nothing
		return sim
    end

	function Base.show(io::IO,sim::Simulation)
		if sim.smoothed[]
	        print(io,"Smoothed Simulation with ", length(sim.domains) ," domain")
		else
			print(io,"Unsmoothed Simulation with ", length(sim.domains)," domain")
		end
		length(sim.domains)>1 ? println(io,"s") : println(io)
        println(io,"\tsites: ",sim.nsites)
        println(io,"\tinterior: ", sim.n)
		println(io,"\texterior: ",sum(sim.exterior))
		println(io,"\tbulk: ",sum(sim.bulk))
        println(io,"\tsurface: ",sum(sim.surface))
		println(io,"\tcorner: ",sum(sim.corner))
		println(io,"\tlattice wall: ",sum(sim.lattice_wall))
		println(io,"\tdomain wall: ",sum(sim.domain_wall))
		println(io,"\t\t====================")
		for d ∈ eachindex(sim.domains)
			println(io,"\t\tDomain ", d,": ",sim.domains[d].type)
		end
	end
end


"""
	smooth!(sim,[num_sub_pixel])

Use sub-pixel sampling along domain walls to smooth the dieletric and pump profiles.

`num_sub_pixel` defaults to `NUM_SUB_PIXEL` in module `Defaults`, which is currently set to 5.
"""
function smooth!(sim::Simulation,num_sub_pixel::Int=NUM_SUB_PIXEL)
    domains = sim.domains
    domain_wall = sim.domain_wall
	lattices = map(z->z.lattice,sim.domains)
    X = Array{Float64}(undef,num_sub_pixel,num_sub_pixel)
    Y = Array{Float64}(undef,num_sub_pixel,num_sub_pixel)
    E = Array{ComplexF64}(undef,num_sub_pixel,num_sub_pixel)
    F = Array{Float64}(undef,num_sub_pixel,num_sub_pixel)
	dw_inds = findall(domain_wall)
    for dw ∈ eachindex(dw_inds)
		d = dw_inds[dw]
		ds = sim.domain_index[d]
	    I,J = Tuple(sim.ij[d])
        is = LinRange(I-.5,I+.5,num_sub_pixel)
        js = LinRange(J-.5,J+.5,num_sub_pixel)
        for i ∈ eachindex(is), j ∈ eachindex(js)
            X[i,j],Y[i,j] = lattices[ds](is[i],js[j])
        end
        simulation_dielectric!(E,domains,X,Y)
        simulation_pump!(F,domains,X,Y)
        sim.ε[d] = mean(E)
        sim.F[d] = mean(F)
    end
	sim.smoothed[] = true
    return nothing
end


# take individually defined domains and concatenate all relevant data
function append_domains!(domain_index,nnxm,nnxp,nnym,nnyp,name,x,y,r,θ,ij,ε,F,interior,bulk,surface,corner,domains)
	count = 1
	for d ∈ domains
		append!(domain_index,fill(count,length(d.x)))
		append!(nnxm,d.nnxm.+length(nnxm))
		append!(nnxp,d.nnxp.+length(nnxp))
		append!(nnym,d.nnym.+length(nnym))
		append!(nnyp,d.nnyp.+length(nnyp))
		append!(name,fill(d.name,length(d.x)))
		append!(x,d.x)
		append!(y,d.y)
		append!(r,d.r)
		append!(θ,d.θ)
		append!(ij,d.ij)
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
get_sites_to_remove(domains,dom,x,y) = which_domains(domains,x,y).!==dom


# identify which domain is associated with a given spatial position (x,y)
function which_domains(domains,x,y)
    doms = Array{Int}(undef,length(x))
	shapes = map(z->z.boundary.shape,domains)
    which_domains!(doms,domains,shapes,x,y)
    return doms
end
function which_domains!(doms::Array{Int},domains,shapes,x,y)
    for i ∈ eachindex(x)
		doms[i] = which_domain(domains,shapes,x[i],y[i])
    end
    return nothing
end
function which_domain(domains,x,y)
	bnd = map(z->z.boundary,domains)
	shapes = map(z->z.boundary.shape,domains)
    return which_domain(domains,shapes,x,y)
end
@inline function which_domain(domains,shapes,x,y)
	dom = findfirst(map(z->z(x,y),shapes))
    isnothing(dom) ? dom=0 : nothing
    iszero(dom) ? nothing : (isvoid(domains[dom]) ? dom=0 : nothing)
    return dom
end
isvoid(domain::Domain) = domain.type==:Void


# remove the pre-computed nearest neighbors, which makes the bulk nn calculation very fast
function trim_nn!(nnxm,nnxp,nnym,nnyp,removed,bulk)
	cr = cumsum(removed)
	for i ∈ eachindex(nnxm)
		if bulk[i]
			nnxm[i] -= cr[nnxm[i]]
			nnxp[i] -= cr[nnxp[i]]
			nnym[i] -= cr[nnym[i]]
			nnyp[i] -= cr[nnyp[i]]
		end
	end
end


# remove the sites determined by get_sites_to_remove
remove_sites!(sites,removed) = map(z->deleteat!(z,findall(removed)),sites)


# compute which sites form the walls along a lattice mismatch and domain walls
function generate_lattice_and_domain_wall!(interior,bulk,surface,corner,domains,domain_index,x,y,ij)
	lattices = map(z->z.lattice,domains)
	shapes = map(z->z.boundary.shape,domains)
	# DX, DY = map(z->z.dx,lattices), map(z->z.dy,lattices)
	domain_wall = falses(size(x)...)
	lattice_wall = falses(size(x)...)
	for i ∈ eachindex(x)
		if interior[i]
			ds = domain_index[i]
			xt,yt = lattices[ds](ij[i]-CartesianIndex(1,0))
			dom1 = which_domain(domains,shapes,xt,yt)
			xt,yt = lattices[ds](ij[i]+CartesianIndex(1,0))
			dom2 = which_domain(domains,shapes,xt,yt)
			xt,yt = lattices[ds](ij[i]-CartesianIndex(0,1))
			dom3 = which_domain(domains,shapes,xt,yt)
			xt,yt = lattices[ds](ij[i]+CartesianIndex(0,1))
			dom4 = which_domain(domains,shapes,xt,yt)

			domain_wall[i] = !(domain_index[i]==dom1==dom2==dom3==dom4)
			domain_wall[i] = domain_wall[i] && (dom1!==0 && dom2!==0 && dom3!==0 && dom4!==0)
			if domain_wall[i]
				lat1 = lattices[dom1]
				lat2 = lattices[dom2]
				lat3 = lattices[dom3]
				lat4 = lattices[dom4]
				lattice_wall[i] = !all(isequal(lattices[domain_index[i]]),(lat1,lat2,lat3,lat4))
				corner[i] = false
				surface[i] = false
				bulk[i] = true
			end
			if lattice_wall[i]
				corner[i] = false
				surface[i] = false
				bulk[i] = false
			end
		end
	end
	return lattice_wall, domain_wall
end


# do a mirror reflection of the surface about the boundary, used in stabilizing the computation of nearest neighbors via a tessellation
function generate_exterior!(corner,domains,domain_index,x,y,surface,verbose::Bool)
	bnds = map(z->z.boundary,domains) ; bcs = map(z->z.bcs,bnds)
	shapes = map(z->z.shape,bnds)
	dxs, dys = map(z->z.lattice.dx,domains), map(z->z.lattice.dy,domains)

	ex_x = NaN*Array{Float64}(undef,sum(surface)+3sum(corner))
	ex_y = NaN*Array{Float64}(undef,sum(surface)+3sum(corner))
	ex_from = zeros(Int,sum(surface)+3sum(corner))
	ex_side = zeros(Int,sum(surface)+3sum(corner))
	normal = Array{Array{Array{Float64,1},1}}(undef,length(x))
	tangent = Array{Array{Array{Float64,1},1}}(undef,length(x))
	distance = Array{Array{Float64,1}}(undef,length(x))
	side = Array{Array{Int,1}}(undef,length(x))

	pg = Progress(sum(surface),.1,"\t\t")
	pg_ref = Threads.Atomic{Int}(0)
	surface_inds = findall(surface .& .!corner)
	for si ∈ eachindex(surface_inds)
		i = surface_inds[si]
		ds = domain_index[i]
        normal[i],tangent[i],distance[i],side[i] = normal_distance(shapes[ds],x[i],y[i])
		ex_x[si] = x[i] - normal[i][1][1]*2distance[i][1]
		ex_y[si] = y[i] - normal[i][1][2]*2distance[i][1]
		ex_from[si] = i
		ex_side[si] = side[i][1]
		verbose ? Threads.atomic_add!(pg_ref,1) : nothing
		if verbose
			update!(pg,pg_ref[])
		end
	end
	corner_inds = findall(corner)
	count = length(surface_inds)+1
	for ci ∈ eachindex(corner_inds)
		i = corner_inds[ci]
		ds = domain_index[i]
        normal[i],tangent[i],distance[i],side[i] = normal_distance(shapes[ds],x[i],y[i])
		dx = min(dxs[ds],dys[ds])
        ex_x[count] = x[i] - normal[i][1][1]*2distance[i][1]
		ex_y[count] = y[i] - normal[i][1][2]*2distance[i][1]
		ex_from[count] = i
		ex_side[count] = side[i][1]
		count += 1
		if corner[i]
			ex_x[count] = x[i] - normal[i][2][1]*2distance[i][2]
			ex_y[count] = y[i] - normal[i][2][2]*2distance[i][2]
			ex_from[count] = i
			ex_side[count] = side[i][2]
			count += 1

			crn_distance = Inf
			crn_index = 0
			for j ∈ eachindex(shapes[ds].corners)
				crn = shapes[ds].corners[j]
				crn_distance_temp = hypot(x[i]-crn[1],y[i]-crn[2])
				if crn_distance_temp < crn_distance
					crn_index=j
					crn_distance = crn_distance_temp
				end
			end
			crn_normal = [x[i]-shapes[ds].corners[crn_index][1],y[i]-shapes[ds].corners[crn_index][2]]
			if norm(crn_normal) > eps()
				crn_normal = crn_normal/norm(crn_normal)
				ex_x[count] = x[i] - crn_normal[1]*2crn_distance
				ex_y[count] = y[i] - crn_normal[2]*2crn_distance
				ex_from[count] = i
				ex_side[count] = length(shapes[ds].corners)+crn_index
				count += 1
			end
		end
	end
	verbose ? finish!(pg) : nothing
	return ex_x[1:count-1], ex_y[1:count-1], ex_from[1:count-1], normal, tangent, distance, side, ex_side[1:count-1]
end


# add the exterior points to the data previously concatenated from the given domains
function append_exterior!(domain_index,nnxm,nnxp,nnym,nnyp,name,x,y,r,θ,ij,ε,F,interior,bulk,surface,corner,
		exterior,domain_wall,lattice_wall,normal,tangent,distance,side,exterior_from,
		exterior_side,exterior_x,exterior_y)

	append!(domain_index,zeros(Int,length(exterior_x)))
	append!(nnxm,zeros(Int,length(exterior_x)))
	append!(nnxp,zeros(Int,length(exterior_x)))
	append!(nnym,zeros(Int,length(exterior_x)))
	append!(nnyp,zeros(Int,length(exterior_x)))
	append!(name,Array{Symbol}(undef,length(exterior_x)))
	append!(ij,Array{CartesianIndex{2}}(undef,length(exterior_x)))
	append!(ε,NaN*ones(ComplexF64,length(exterior_x)))
	append!(F,zeros(Float64,length(exterior_x)))

	append!(interior,falses(length(exterior_x)))
	append!(bulk,falses(length(exterior_x)))
	append!(surface,falses(length(exterior_x)))
	append!(corner,falses(length(exterior_x)))
	append!(exterior,trues(length(exterior_x)))
	append!(domain_wall,falses(length(exterior_x)))
	append!(lattice_wall,falses(length(exterior_x)))

	append!(normal,Array{Array{Array{Float64,1},1}}(undef,length(exterior_x)))
	append!(tangent,Array{Array{Array{Float64,1},1}}(undef,length(exterior_x)))
	append!(distance,Array{Array{Float64,1}}(undef,length(exterior_x)))
	append!(side,Array{Array{Int,1}}(undef,length(exterior_x)))

	prepend!(exterior_from,zeros(Int,length(x)))
	prepend!(exterior_side,zeros(Int,length(x)))
	append!(x,exterior_x)
	append!(y,exterior_y)
	append!(r,NaN*ones(Float64,length(exterior_x)))
	return nothing
end


# KEY BLOCKS: LINKS BULK POINTS via precomputed nn and lookup tables
function link_bulk!(link_fiduciary::BitArray,sim::Simulation,verbose::Bool)

	domains,x,y,ij = sim.domains, sim.x, sim.y, sim.ij
	lattices = map(z->z.lattice,domains)

	DX, DY = map(z->z.dx,lattices), map(z->z.dy,lattices)
	V1, V2 = map(z->z.v1,lattices), map(z->z.v2,lattices)
	DR, DΘ = map(z->z.dr,lattices), map(z->z.dθ,lattices)

	ijlat = map((x,y)->(x,lattices[y]),ij[sim.interior],sim.domain_index[sim.interior])

	pg = Progress(sum(sim.bulk),.1,"\t\t")
	pg_ref = Threads.Atomic{Int}(0)
	c0 = 0
	bulk_inds = findall(sim.bulk)
	# Threads.@threads
	for bi ∈ eachindex(bulk_inds)
		c = c0 + 4(bi-1)+1
		i = bulk_inds[bi]
        begin
			ds = sim.domain_index[i]
			lds = lattices[ds]
			dx, dy = DX[ds], DY[ds]
			dr, dθ = DR[ds], DΘ[ds]

			sim.link_start[c] = i
			sim.link_start[c+1] = i
			sim.link_start[c+2] = i
			sim.link_start[c+3] = i

			if !sim.domain_wall[i]
				sim.link_stop[c  ] = sim.nnxm[i]
	            sim.link_stop[c+1] = sim.nnxp[i]
	            sim.link_stop[c+2] = sim.nnym[i]
	            sim.link_stop[c+3] = sim.nnyp[i]
			else
				sim.link_stop[c  ] = findfirst(isequal((ij[i]+CartesianIndex(-1, 0),lds)),ijlat)
				sim.link_stop[c+1] = findfirst(isequal((ij[i]+CartesianIndex( 1, 0),lds)),ijlat)
				sim.link_stop[c+2] = findfirst(isequal((ij[i]+CartesianIndex( 0,-1),lds)),ijlat)
				sim.link_stop[c+3] = findfirst(isequal((ij[i]+CartesianIndex( 0, 1),lds)),ijlat)
			end

			if lds.type==:Cartesian
				sim.link_weight[c  ] = 1/dx^2
				sim.link_weight[c+1] = 1/dx^2
				sim.link_weight[c+2] = 1/dy^2
				sim.link_weight[c+3] = 1/dy^2
			elseif lds.type==:Polar
				sim.link_weight[c  ] = 1/dr^2
				sim.link_weight[c+1] = 1/dr^2
				sim.link_weight[c+2] = 1/dθ^2
				sim.link_weight[c+3] = 1/dθ^2

				sim.link_half_r[c] = lds.r0 + (sim.ij[i][1]-1/2)*lds.dr
				sim.link_half_r[c+1] = lds.r0 + (sim.ij[i][1]+1/2)*lds.dr
				sim.link_r_bool[c] = true
			else
				throw(LatticeError(lds.type))
			end

			sim.link_half_x[c],sim.link_half_y[c] = lds(Tuple(sim.ij[i]).-(1/2,0))
			sim.link_x_bool[c] = true

			sim.link_half_x[c+1],sim.link_half_y[c+1] = lds(Tuple(sim.ij[i]).+(1/2,0))

			sim.link_x_bool[c+1] = true

			sim.link_half_x[c+2],sim.link_half_y[c+2] = lds(Tuple(sim.ij[i]).-(0,1/2))
			sim.link_y_bool[c+2] = true

			sim.link_half_x[c+3],sim.link_half_y[c+3] = lds(Tuple(sim.ij[i]).+(0,1/2))
			sim.link_y_bool[c+3] = true

			link_fiduciary[c] = true
			link_fiduciary[c+1] = true
			link_fiduciary[c+2] = true
			link_fiduciary[c+3] = true

		end
		verbose ? Threads.atomic_add!(pg_ref,1) : nothing
		if verbose && Threads.threadid()==Threads.nthreads()
			update!(pg,pg_ref[])
		end
    end
	verbose ? finish!(pg) : nothing
    return nothing
end


# KEY BLOCKS: LINKS LATTICE WALL POINTS via interpolation
function link_lattice_wall!(link_fiduciary::BitArray,sim::Simulation,verbose::Bool)

	NNN = NNN_BULK

	domains,x,y,ij = sim.domains,sim.x,sim.y,sim.ij

	lattices = map(z->z.lattice,domains)
	DX, DY = map(z->z.dx,lattices), map(z->z.dy,lattices)
	DR, DΘ = map(z->z.dr,lattices), map(z->z.dθ,lattices)

	pg = Progress(sum(sim.lattice_wall),.1,"\t\t")
	pg_ref = Threads.Atomic{Int}(0)
	c0 = findlast(link_fiduciary)
	isnothing(c0) ? c0=0 : nothing

	lattice_wall_inds = findall(sim.lattice_wall)
	# Threads.@threads
	for lwi ∈ eachindex(lattice_wall_inds)
		c = c0 + 4NNN*(lwi-1)+1
		i = lattice_wall_inds[lwi]
		begin
			ds = sim.domain_index[i]
			lds = lattices[ds]

			if lds.type==:Cartesian
				dx, dy = DX[ds], DY[ds]
			elseif lds.type==:Polar
				dr, dθ = DR[ds], DΘ[ds]
				dx = dr
				dy = sim.r[i]*dθ
			else
				throw(LatticeError(lds.type))
			end

			x1, y1 = lds(Tuple(sim.ij[i]).-(1,0))
			x2, y2 = lds(Tuple(sim.ij[i]).+(1,0))
			x3, y3 = lds(Tuple(sim.ij[i]).-(0,1))
			x4, y4 = lds(Tuple(sim.ij[i]).+(0,1))

			hx1, hy1 = lds(Tuple(sim.ij[i]).-(1/2,0))
			hx2, hy2 = lds(Tuple(sim.ij[i]).+(1/2,0))
			hx3, hy3 = lds(Tuple(sim.ij[i]).-(0,1/2))
			hx4, hy4 = lds(Tuple(sim.ij[i]).+(0,1/2))

			it,wt = generate_inds_weights(x1,y1,sim,NNN,i,ds)
			N = length(it)
			for j ∈ 1:N
				sim.link_start[c+j-1] = i
				sim.link_stop[c+j-1] = it[j]
				sim.link_weight[c+j-1] = wt[j]/dx^2
				sim.link_half_x[c+j-1] = hx1
				sim.link_half_y[c+j-1] = hy1
				sim.link_x_bool[c+j-1] = true
				link_fiduciary[c+j-1] = true
			end
			c += NNN

			it,wt = generate_inds_weights(x2,y2,sim,NNN,i,ds)
			N = length(it)
			for j ∈ 1:N
				sim.link_start[c+j-1] = i
				sim.link_stop[c+j-1] = it[j]
				sim.link_weight[c+j-1] = wt[j]/dx^2
				sim.link_half_x[c+j-1] = hx2
				sim.link_half_y[c+j-1] = hy2
				sim.link_x_bool[c+j-1] = true
				link_fiduciary[c+j-1] = true
			end
			c += NNN

			it,wt = generate_inds_weights(x3,y3,sim,NNN,i,ds)
			N = length(it)
			for j ∈ 1:N
				sim.link_start[c+j-1] = i
				sim.link_stop[c+j-1] = it[j]
				sim.link_weight[c+j-1] = wt[j]/dy^2
				sim.link_half_x[c+j-1] = hx3
				sim.link_half_y[c+j-1] = hy3
				sim.link_y_bool[c+j-1] = true
				link_fiduciary[c+j-1] = true
			end
			c += NNN

			it,wt = generate_inds_weights(x4,y4,sim,NNN,i,ds)
			N = length(it)
			for j ∈ 1:N
				sim.link_start[c+j-1] = i
				sim.link_stop[c+j-1] = it[j]
				sim.link_weight[c+j-1] = wt[j]/dy^2
				sim.link_half_x[c+j-1] = hx4
				sim.link_half_y[c+j-1] = hy4
				sim.link_y_bool[c+j-1] = true
				link_fiduciary[c+j-1] = true
			end
		end
		verbose ? Threads.atomic_add!(pg_ref,1) : nothing
		if verbose && Threads.threadid()==Threads.nthreads()
			update!(pg,pg_ref[])
		end
	end
	verbose ? finish!(pg) : nothing
	return nothing
end


# KEY BLOCKS: LINKS SURFACE POINTS via interpolation
function link_surface!(link_fiduciary::BitArray,self_fiduciary::BitArray,sim::Simulation,verbose::Bool)

	domains,x,y = sim.domains,sim.x,sim.y
	lattices = map(z->z.lattice,domains)
	DX, DY = map(z->z.dx,lattices), map(z->z.dy,lattices)
	DR, DΘ = map(z->z.dr,lattices), map(z->z.dθ,lattices)

	surface_inds = findall(sim.surface .& .!sim.corner)
	pg = Progress(length(surface_inds),.1,"\t\t")
	pg_ref = Threads.Atomic{Int}(0)
	cl0 = findlast(link_fiduciary)
	isnothing(cl0) ? cl0=0 : nothing
	cs0 = 0
	# Threads.@threads
	for si ∈ eachindex(surface_inds)

		i = surface_inds[si]
		cs = cs0 + si
		cl = cl0 + (NNN_BULK+2NNN_SURFACE)*(si-1) + 1

		begin
			ds = sim.domain_index[i]
			n,t,d,sides = sim.normal[i],sim.tangent[i],sim.distance[i],sim.side[i]
			lds = lattices[ds]

			dxn = 2d[1]
			NNN = NNN_BULK

			hxn1, hyn1 = x[i] - n[1][1]*dxn/2, y[i] - n[1][2]*dxn/2
			sim.self_index[cs] = i
			sim.self_weight[cs] = -1/dxn^2
			sim.self_half_x[cs] = hxn1
			sim.self_half_y[cs] = hyn1
			if sides[1] ∈ (1,2)
				sim.self_x_bool[cs] = true
			else
				sim.self_y_bool[cs] = true
			end
			self_fiduciary[cs] = true

			if lds.type==:Cartesian
				dxn = max(dxn,min(DX[ds],DY[ds]))
			elseif lds.type==:Polar
				dxn = max(dxn,min(DR[ds],sim.r[i]*DΘ[ds]))
			else
				throw(LatticeError(lds.type))
			end

			xn2, yn2   = x[i] + n[1][1]*dxn  , y[i] + n[1][2]*dxn
			hxn2, hyn2 = x[i] + n[1][1]*dxn/2, y[i] + n[1][2]*dxn/2
			it,wt = generate_inds_weights(xn2,yn2,sim,NNN,i,ds)
			N = length(it)
			for j ∈ 1:N
				sim.link_start[cl+j-1] = i
				sim.link_stop[cl+j-1] = it[j]
				sim.link_weight[cl+j-1] = wt[j]/dxn^2
				sim.link_half_x[cl+j-1] = hxn2
				sim.link_half_y[cl+j-1] = hyn2
				if sides[1] ∈ (1,2)
					sim.link_x_bool[cl+j-1] = true
				else
					sim.link_y_bool[cl+j-1] = true
				end
				link_fiduciary[cl+j-1] = true
			end
			cl += NNN

			if lds.type==:Cartesian
				dxt = min(DX[ds],DY[ds])
			elseif lds.type==:Polar
				dxt = min(DR[ds],sim.r[i]*DΘ[ds])
			else
				throw(LatticeError(lds.type))
			end

			NNN =  NNN_SURFACE

			xt1, yt1   = x[i] - t[1][1]*dxt  , y[i] - t[1][2]*dxt
			hxt1, hyt1 = x[i] - t[1][1]*dxt/2, y[i] - t[1][2]*dxt/2

			it,wt = generate_inds_weights(xt1,yt1,sim,NNN,i,ds)
			N = length(it)
			for j ∈ 1:N
				sim.link_start[cl+j-1] = i
				sim.link_stop[cl+j-1] = it[j]
				sim.link_weight[cl+j-1] = wt[j]/dxt^2
				sim.link_half_x[cl+j-1] = hxt1
				sim.link_half_y[cl+j-1] = hyt1
				if sides[1] ∈ (1,2)
					sim.link_y_bool[cl+j-1] = true
				else
					sim.link_x_bool[cl+j-1] = true
				end
				link_fiduciary[cl+j-1] = true
			end
	        cl += NNN

			xt2, yt2   = x[i] + t[1][1]*dxt  , y[i] + t[1][2]*dxt
			hxt2, hyt2 = x[i] + t[1][1]*dxt/2, y[i] + t[1][2]*dxt/2
			it,wt = generate_inds_weights(xt2,yt2,sim,NNN,i,ds)
			N = length(it)
			for j ∈ 1:N
				sim.link_start[cl+j-1] = i
				sim.link_stop[cl+j-1] = it[j]
				sim.link_weight[cl+j-1] = wt[j]/dxt^2
				sim.link_half_x[cl+j-1] = hxt2
				sim.link_half_y[cl+j-1] = hyt2
				if sides[1] ∈ (1,2)
					sim.link_y_bool[cl+j-1] = true
				else
					sim.link_x_bool[cl+j-1] = true
				end
				link_fiduciary[cl+j-1] = true
			end
		end
		verbose ? Threads.atomic_add!(pg_ref,1) : nothing
		if verbose && Threads.threadid()==Threads.nthreads()
			update!(pg,pg_ref[])
		end
	end
	verbose ? finish!(pg) : nothing
    return nothing
end


# KEY BLOCKS: LINKS CORNER POINTS via interpolation
function link_corner!(link_fiduciary::BitArray,self_fiduciary::BitArray,sim::Simulation,verbose::Bool)
	domains,x,y = sim.domains,sim.x,sim.y

	NNN = NNN_CORNER

	cl0 = findlast(link_fiduciary)
	cs0 = findlast(self_fiduciary)

	lattices = map(z->z.lattice,domains)
	DX, DY = map(z->z.dx,lattices), map(z->z.dy,lattices)
	DR, DΘ = map(z->z.dr,lattices), map(z->z.dθ,lattices)

	corner_inds = findall(sim.corner)
	for ci ∈ eachindex(corner_inds)
		i = corner_inds[ci]
		cl = cl0 + 2NNN*(ci-1) + 1
		cs = cs0 + 2(ci-1) + 1
		begin
			ds = sim.domain_index[i]
			n,t,d,sides = sim.normal[i],sim.tangent[i],sim.distance[i],sim.side[i]
			lds = lattices[ds]

			dxn = 2d[1]

			hxn1, hyn1 = x[i] - n[1][1]*dxn/2, y[i] - n[1][2]*dxn/2
			sim.self_index[cs] = i
			sim.self_weight[cs] = -1/dxn^2
			sim.self_half_x[cs] = hxn1
			sim.self_half_y[cs] = hyn1
			self_fiduciary[cs] = true
			if sides[1] ∈ (1,2)
				sim.self_x_bool[cs] = true
			else
				sim.self_y_bool[cs] = true
			end
			cs += 1


			if lds.type==:Cartesian
				dxn = max(dxn,min(DX[ds],DY[ds]))
			elseif lds.type==:Polar
				dxn = max(dxn,min(DR[ds],DΘ[ds]))
			else
				throw(LatticeError(lds.type))
			end

			xn2, yn2   = x[i] + n[1][1]*dxn  , y[i] + n[1][2]*dxn
			hxn2, hyn2 = x[i] + n[1][1]*dxn/2, y[i] + n[1][2]*dxn/2
			it,wt = generate_inds_weights(xn2,yn2,sim,NNN,i,ds)
			N = length(it)
			for j ∈ 1:N
				sim.link_start[cl+j-1] = i
				sim.link_stop[cl+j-1] = it[j]
				sim.link_weight[cl+j-1] = wt[j]/dxn^2
				sim.link_half_x[cl+j-1] = hxn2
				sim.link_half_y[cl+j-1] = hyn2
				link_fiduciary[cl+j-1] = true
				if sides[1] ∈ (1,2)
					sim.link_x_bool[cl+j-1] = true
				else
					sim.link_y_bool[cl+j-1] = true
				end
			end
			cl += NNN

			dxn = 2d[2]

			hxn1, hyn1 = x[i] - n[2][1]*dxn/2, y[i] - n[2][2]*dxn/2
			sim.self_index[cs] = i
			sim.self_weight[cs] = -1/dxn^2
			sim.self_half_x[cs] = hxn1
			sim.self_half_y[cs] = hyn1
			self_fiduciary[cs] = true
			if sides[1] ∈ (1,2)
				sim.self_y_bool[cs] = true
			else
				sim.self_x_bool[cs] = true
			end

			if lds.type==:Cartesian
				dxn = max(dxn,min(DX[ds],DY[ds]))
			elseif lds.type==:Polar
				dxn = max(dxn,min(DR[ds],DΘ[ds]))
			else
				throw(LatticeError(lds.type))
			end

			xn2, yn2   = x[i] + n[2][1]*dxn  , y[i] + n[2][2]*dxn
			hxn2, hyn2 = x[i] + n[2][1]*dxn/2, y[i] + n[2][2]*dxn/2
			it,wt = generate_inds_weights(xn2,yn2,sim,NNN,i,ds)
			N = length(it)
			for j ∈ 1:N
				sim.link_start[cl+j-1] = i
				sim.link_stop[cl+j-1] = it[j]
				sim.link_weight[cl+j-1] = wt[j]/dxn^2
				sim.link_half_x[cl+j-1] = hxn2
				sim.link_half_y[cl+j-1] = hyn2
				link_fiduciary[cl+j-1] = true
				if sides[1] ∈ (1,2)
					sim.link_y_bool[cl+j-1] = true
				else
					sim.link_x_bool[cl+j-1] = true
				end
			end
		end
	end
    return nothing
end


# trim the preallocated link arrays to the size actually used
function trim_link!(sim::Simulation,link_fiduciary::BitArray)
	deleteat!(sim.link_start,.!link_fiduciary)
	deleteat!(sim.link_stop,.!link_fiduciary)
	deleteat!(sim.link_weight,.!link_fiduciary)
	deleteat!(sim.link_half_x,.!link_fiduciary)
	deleteat!(sim.link_half_y,.!link_fiduciary)
	deleteat!(sim.link_half_r,.!link_fiduciary)

	link_x_bool = Array(sim.link_x_bool)
	deleteat!(link_x_bool,.!link_fiduciary)
	sim.link_x_bool[1:length(link_x_bool)] = link_x_bool
	resize!(sim.link_x_bool,length(link_x_bool))

	link_y_bool = Array(sim.link_y_bool)
	deleteat!(link_y_bool,.!link_fiduciary)
	sim.link_y_bool[1:length(link_y_bool)] = link_y_bool
	resize!(sim.link_y_bool,length(link_y_bool))

	link_r_bool = Array(sim.link_r_bool)
	deleteat!(link_r_bool,.!link_fiduciary)
	sim.link_r_bool[1:length(link_r_bool)] = link_r_bool
	resize!(sim.link_r_bool,length(link_r_bool))

	return nothing
end


# trim the preallocated self-linking arrays to the size actually used
function trim_self!(sim::Simulation,self_fiduciary::BitArray)
	deleteat!(sim.self_index,.!self_fiduciary)
	deleteat!(sim.self_weight,.!self_fiduciary)
	deleteat!(sim.self_half_x,.!self_fiduciary)
	deleteat!(sim.self_half_y,.!self_fiduciary)
	deleteat!(sim.self_half_r,.!self_fiduciary)

	self_x_bool = Array(sim.self_x_bool)
	deleteat!(self_x_bool,.!self_fiduciary)
	sim.self_x_bool[1:length(self_x_bool)] = self_x_bool
	resize!(sim.self_x_bool,length(self_x_bool))

	self_y_bool = Array(sim.self_y_bool)
	deleteat!(self_y_bool,.!self_fiduciary)
	sim.self_y_bool[1:length(self_y_bool)] = self_y_bool
	resize!(sim.self_y_bool,length(self_y_bool))

	self_r_bool = Array(sim.self_r_bool)
	deleteat!(self_r_bool,.!self_fiduciary)
	sim.self_r_bool[1:length(self_r_bool)] = self_r_bool
	resize!(sim.self_r_bool,length(self_r_bool))

	return nothing
end


# apply the local boundary conditions (Dirichlet, Neumann, etc)
function apply_local_bc!(sim::Simulation)
	apply_local_bc_link!(sim)
	apply_local_bc_self!(sim)
end

# apply local boundary conditions to the links
function apply_local_bc_link!(sim::Simulation)
	for i ∈ eachindex(sim.link_start)
		if sim.exterior[sim.link_stop[i]]
			new_stop = sim.exterior_from[sim.link_stop[i]]
			side = sim.exterior_side[sim.link_stop[i]]
			bcs = sim.domains[sim.domain_index[new_stop]].boundary.bcs
			shp = sim.domains[sim.domain_index[new_stop]].boundary.shape
			if typeof(shp)<:AbstractQuadrilateral
				if side≤length(shp.corners)
					if typeof(bcs[side])<:DirichletBC
						sim.link_weight[i] *= -1
		            elseif typeof(bcs[side])<:NeumannBC
		                nothing
					elseif typeof(bcs[side])<:noBC
						sim.link_weight[i] *= 0
		            end
				else
					if side-length(shp.corners)==1
						sides = [1,3]
					elseif side-length(shp.corners)==2
						sides = [2,3]
					elseif side-length(shp.corners)==3
						sides = [1,4]
					else
						sides = [2,4]
					end
					for j ∈ sides
						if typeof(bcs[j])<:DirichletBC
							sim.link_weight[i] *= -1
						elseif typeof(bcs[j])<:NeumannBC
							nothing
						elseif typeof(bcs[j])<:noBC
							sim.link_weight[i] *= 0
						end
					end
				end
			else
				if typeof(bcs[side])<:DirichletBC
					sim.link_weight[i] *= -1
	            elseif typeof(bcs[side])<:NeumannBC
	                nothing
				elseif typeof(bcs[side])<:noBC
					sim.link_weight[i] *= 0
	            end
			end
			sim.link_stop[i] = new_stop
		end
	end
end

# apply local boundary conditions to the self-links
function apply_local_bc_self!(sim::Simulation)
	for i ∈ eachindex(sim.self_index)
		si = sim.self_index[i]
		if !sim.corner[si]
			side = sim.side[si][1]
			bcs = sim.domains[sim.domain_index[si]].boundary.bcs
			if typeof(bcs[side])<:DirichletBC
				sim.self_weight[i] *= 2
	        elseif typeof(bcs[side])<:NeumannBC
	            sim.self_weight[i] *= 0
			elseif typeof(bcs[side])<:noBC
				nothing
	        end
		else
			if sim.self_x_bool[i]
				side = sim.side[si][1]
			else
				side = sim.side[si][2]
			end
			bcs = sim.domains[sim.domain_index[si]].boundary.bcs
			if typeof(bcs[side])<:DirichletBC
				sim.self_weight[i] *= 2
			elseif typeof(bcs[side])<:NeumannBC
				sim.self_weight[i] *= 0
			elseif typeof(bcs[side])<:noBC
				nothing
			end
		end
	end
end


# apply the nonlocal boundary conditions (not implemented yet)
function apply_nonlocal_bc!(sim::Simulation)
	return nothing
end


# get the interpolation indices and weights for (X,Y) on sites (x,y), with a maximal interpolation order of `NNN`
# do this by first finding the surrounding sites 2 layers deep (from Tessellations) and then finding the subset
# of N points which give the best-conditioned matrix and computing the weights for those points.
@inline function generate_inds_weights(X::Real,Y::Real,sim::Simulation,NNN::Int,site_index::Int,domain_index::Int)
	# t = sim.tessellation
	NNT = ceil(Int,sqrt(N_NEAREST_SITES_FOR_TESSELLATION/π))
	dx = sim.domains[domain_index].lattice.dx
	dy = sim.domains[domain_index].lattice.dy
	inds = hypot.(sim.x.-X,sim.y.-Y) .< max(dx,dy)*NNT
	deep = false
	tess_inds = get_surrounding_sites(Tessellation(sim.x[inds],sim.y[inds]),X,Y,deep)
	inds = findall(inds)[tess_inds]

	# inds = inds[hypot.(X.-x[inds],Y.-y[inds]).≤2NNN*max(dx,dy)]
	# inds = inds[sim.interior[inds]]
	# surface && (sum(map(z->z.boundary.shape(x[site_index],y[site_index]),sim.domains)) ≤ 1) ? setdiff!(inds,sim.not_this_domains_exterior[domain_index]) : nothing
	# surface=true
	# surface ?
	# setdiff!(inds,sim.not_this_domains_exterior[domain_index])# : setdiff!(inds,sim.exterior_inds)

	NN = min(length(inds),NNN)
	if N_CUBIC ≤ NN
		N = N_CUBIC
	elseif N_QUADRATIC ≤ NN
		N = N_QUADRATIC
	elseif N_BILINEAR ≤ NN
		N = N_BILINEAR
	elseif N_LINEAR ≤ NN
		N = N_LINEAR
	end
	# perm = sortperm(hypot.(X.-sim.x[inds],Y.-sim.y[inds]))
	# inds = inds[perm[1:N]]

	# find the combination of indices with the smallest condition number
	x = sim.x; y=sim.y
	C = collect(combinations(inds,N))
	conds = Ref(Inf); ii = Ref(1)
	for i ∈ eachindex(C)
		temp = cond(generate_interpolating_matrix(X,Y,x[C[i]],y[C[i]]))
		if temp < conds[]
			conds[] = temp
			ii[] = i
		end
	end
	inds = C[ii[]]

    weights = generate_weights(X,Y,x[inds],y[inds])
	return inds, weights
end


const POLYVEC3 = [1,0,0]
const POLYVEC4 = [1,0,0,0]
const POLYVEC6 = [1,0,0,0,0,0]
const POLYVEC10 = [1,0,0,0,0,0,0,0,0,0]
# "invert" the interpolating matrix or do a pseudo-inverse if singular
@inline function generate_weights(X::Real,Y::Real,x::Array,y::Array)
	N = length(x)
	interp_matrix = generate_interpolating_matrix(X,Y,x,y)
	try
		if N==3
			weights = interp_matrix\POLYVEC3
		elseif N==4
			weights = interp_matrix\POLYVEC4
		elseif N==6
			weights = interp_matrix\POLYVEC6
		elseif N==10
			weights = interp_matrix\POLYVEC10
		end
		return weights
	catch
		return weights = pinv(interp_matrix)[:,1]
	end
end


# matrix of powers of relative position, whose inverse contains the interpolation coefficients
@inline function generate_interpolating_matrix(X::Real,Y::Real,x::Array,y::Array)
	N = length(x)
	@assert N≤10 "only coded for N≤10"
	interp_matrix = Matrix{Float64}(undef,N,N)
	@inbounds for i ∈ 1:N
		N≥1  ? interp_matrix[1 ,i] = 1 : nothing
		N≥2  ? interp_matrix[2 ,i] = (y[i]-Y) : nothing
		N≥3  ? interp_matrix[3 ,i] = (x[i]-X) : nothing
		N≥4  ? interp_matrix[4 ,i] = (x[i]-X)*(y[i]-Y) : nothing
		N≥5  ? interp_matrix[5 ,i] = (y[i]-Y)^2 : nothing
		N≥6  ? interp_matrix[6 ,i] = (x[i]-X)^2 : nothing
		N≥7  ? interp_matrix[7 ,i] = (x[i]-X)^1*(y[i]-Y)^2 : nothing
		N≥8  ? interp_matrix[8 ,i] = (x[i]-X)^2*(y[i]-Y)^1 : nothing
		N≥9  ? interp_matrix[9 ,i] = (y[i]-Y)^3 : nothing
		N≥10 ? interp_matrix[10,i] = (x[i]-X)^3 : nothing
	end
	return interp_matrix
end


"""
	update_dielectric!(sim,[ε])

update simulation `sim` by recomputing the dielectric function or by giving it an explicit array `ε`
"""
function update_dielectric!(sim::Simulation,ε::Array)
	sim.ε[:] = ε
	sim.smoothed[] = false
	return nothing
end
update_dielectric!(sim::Simulation) = update_dielectric!(sim,simulation_dielectric.(Ref(sim),sim.x,sim.y))
"""
	update_pump!(sim,[F])

update simulation `sim` by recomputing the pump function or by giving it an explicit array `F`
"""
function update_pump!(sim::Simulation,F::Array)
	sim.F[:] = F
	sim.smoothed[] = false
	return nothing
end
update_pump!(sim::Simulation) = update_pump!(sim,simulation_pump.(Ref(sim),sim.x,sim.y))


"""
	update_sim!(sim)

update simulation `sim` by recomputing the dielectric and pump functions
"""
function update_sim!(sim::Simulation)
	update_dielectric!(sim)
	update_pump!(sim)
	return nothing
end


# used in the smoothing routing smooth!, it computes the dielectric as a function of position (x,y)
function simulation_dielectric!(ε,domains,x,y)
    for i ∈ eachindex(x)
        ε[i] = simulation_dielectric(domains,x[i],y[i])
    end
    return nothing
end
function simulation_dielectric(domains,x,y)
    dom = which_domain(domains,x,y)
    ε = dom==0 ? complex(NaN,NaN) : domains[dom].dielectric(x,y)
    return ε
end
simulation_dielectric(sim::Simulation,x,y) = simulation_dielectric(sim.domains,x,y)


# used in the smoothing routing smooth!, it computes the pump as a function of position (x,y)
function simulation_pump!(F,domains,x,y)
    for i ∈ eachindex(x)
        F[i] = simulation_pump(domains,x[i],y[i])
    end
    return nothing
end
function simulation_pump(domains,x,y)
    dom = which_domain(domains,x,y)
    F = dom==0 ? 0.0 : domains[dom].pump(x,y)
    return F
end
simulation_pump(sim::Simulation,x,y) = simulation_pump(sim.domains,x,y)


####################################################################################
# PLOTTING
####################################################################################

# plotting sim, with by argument, for external use
@recipe function f(sim::Simulation; by=nothing)
	if isnothing(by)
		layout --> (1,3)
	else
		layout --> 1
	end
	(sim, by)
end
@recipe function f(sim::Simulation,by::Union{Symbol,Function,Nothing})
	if isnothing(by)
		bys = [:real,:imag,:F]
		for i ∈ eachindex(bys)
			@series begin
				subplot --> i
				markersize --> MARKERSIZE_SCALE/sqrt(length(sim.x))/3
				(sim,bys[i],1)
			end
		end
	else
		markersize --> MARKERSIZE_SCALE/sqrt(length(sim.x))
		(sim,by,1)
	end
end

# elementary Simulation plot, for internal use
@recipe function f(sim::Simulation,by::Union{Symbol,Function},_1::Int)
	aspect_ratio --> 1
	legend --> false
	if by ∈ [:real,real,:Real,:Re,:re]
		z = real(sqrt.(sim.ε[sim.interior]))
		@series begin
			seriestype --> :scatter
			markershape --> MARKERSHAPE
			markerstrokealpha --> 0
	        clims --> (minimum(z), maximum(z))
			color --> :sequential
			marker_z --> z
			(sim.x[sim.interior], sim.y[sim.interior])
		end
		for i ∈ reverse(eachindex(sim.domains))
			@series begin
				title --> "Re n(x)"
				sim.domains[i].boundary
			end
		end
	elseif by ∈ [:imag,imag,:Imag,:Im,:im]
		z = imag(sqrt.(sim.ε[sim.interior]))
		@series begin
			seriestype --> :scatter
			markershape --> MARKERSHAPE
			markerstrokealpha --> 0
			clims --> (-maximum(abs.(z)), maximum(abs.(z)))
			color --> :diverging
			marker_z --> z
			(sim.x[sim.interior], sim.y[sim.interior])
		end
		for i ∈ reverse(eachindex(sim.domains))
			@series begin
				title --> "Im n(x)"
				sim.domains[i].boundary
			end
		end
	elseif by ∈ [:F,:pump,:f]
		z = sim.F[sim.interior]
		@series begin
			seriestype --> :scatter
			markershape --> MARKERSHAPE
			markerstrokealpha --> 0
			clims --> (-maximum(abs.(z)), maximum(abs.(z)))
			color --> :diverging
	        marker_z --> z
			(sim.x[sim.interior], sim.y[sim.interior])
		end
		for i ∈ reverse(eachindex(sim.domains))
			@series begin
				title --> "F(x)"
				sim.domains[i].boundary
			end
		end
	elseif by ∈ [:links, :link, :Links, :Link, :Hopping, :hopping, :adjacency, :Adjacency]
		@series (sim,true)
	else
		throw("unrecognized `by` keyword $by, must be one of :real, :imag, :F")
	end
end

@recipe function f(sim::Simulation,islink::Bool)
	legend --> false
	aspect_ratio --> 1
	S = sparse(sim.link_start,sim.link_stop,sim.link_weight,length(sim.x),length(sim.x),+)
	I,J,W = findnz(S)
	perm = sortperm(W,by=abs)
	I = I[perm]; J = J[perm]; W = abs.(W[perm])
	for i ∈ eachindex(I)
		inds = findall(isequal(i),I)
		W[inds] = W[inds]/mean(W[inds])
	end
	X = [sim.x[I] sim.x[J] NaN*ones(length(I))]
	Y = [sim.y[I] sim.y[J] NaN*ones(length(I))]
	Z = abs.([W W W])
	for i ∈ reverse(eachindex(sim.domains))
		@series begin
			ms --> LINK_SCALE_REDUCTION*MARKERSIZE_SCALE/sqrt(length(sim.x))
			sim.domains[i]
		end
	end
	@series begin
		clim --> (0,2)
		colorbar --> false
		seriestype --> :line
		aspect_ratio --> 1
		line_z --> Z'
		color--> :sequential_r
		(X',Y')
	end
end



end # module

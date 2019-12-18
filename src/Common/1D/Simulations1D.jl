# TODO: fix and test smoothing

function Simulation(
			ω₀::Real,
			lattice_domains::NTuple{K,LatticeDomain{1}},
			nondispersive_domains::NTuple{L,NondispersiveDomain{1}},
			dispersive_domains::NTuple{M,DispersiveDomain{1}};
			k₂₀::Real=0,
			k₃₀::Real=0,
			k₁₀::Real=sqrt(ω₀^2-k₂₀^2-k₃₀^2)
			) where {K,L,M}

	# simple checking for 1D symmetric
	K==1 || throw("in 1D, can only provide one LatticeDomain, but $K provided instead")
	lattice_domain = lattice_domains[1]

	smoothed = Ref(false)
	x = lattice_domain.x

	lattice_domain_indices = fill(1,length(x))
	nondispersive_domain_indices = which_domains(nondispersive_domains,x)
	dispersive_domain_indices = which_domains(dispersive_domains,x)

	ε = Vector{ComplexF64}(undef,length(x))
	if isempty(nondispersive_domains)
		for i ∈ eachindex(ε) ε[i] = lattice_domain.ε end
	else
		dielectrics = map(n -> getfield(n,:dielectric), nondispersive_domains)
		for i ∈ eachindex(ε)
			d = nondispersive_domain_indices[i]
			if iszero(d)
				ε[i] = lattice_domain.ε
			else
				ε[i] = dielectrics[d](x[i])
			end
		end
	end

	F = Vector{Float64}(undef,length(x))
	χ = Vector{AbstractDispersion}(undef,length(x))
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
				F[i] = pumps[d](x[i])
				χ[i] = χs[d]
			end
		end
	end

	indices = lattice_domain.indices
	nnm = lattice_domain.nnm
	nnp = lattice_domain.nnp
	interior = lattice_domain.interior
	surface = lattice_domain.surface
	bulk = lattice_domain.bulk
	corner = lattice_domain.corner

	# checking
	surface[1] || throw("first site should be an endpoint, something is amiss")
	surface[end] || throw("last site should be an endpoint, something is amiss")
	sum(surface)>2 && throw("no point except first and last should be an endpoint, something is amiss")

	x_half = _generate_half_x(lattice_domain.lattice, x, indices)

	σ = (Vector{ComplexF64}(undef,length(x)),)
	σ_half = (Vector{ComplexF64}(undef,length(x)+1),)
	boundary_layer!(σ, lattice_domains, x, lattice_domain_indices)
	boundary_layer!(σ_half, lattice_domains, x_half[1], vcat(lattice_domain_indices, lattice_domain_indices[end]))
	α = (1 .+ 1im*σ[1]/k₁₀,)
	α_half = (1 .+ 1im*σ_half[1]/k₁₀,)

	laplacian = Laplacian{Symmetric}(lattice_domain.lattice,α[1],α_half[1])
	curlcurl = Curlcurl(lattice_domain.lattice,α[1],α_half[1],nnm,nnp,indices,interior,surface)

	Σ = SelfEnergy(lattice_domain,α_half[1])

	return Simulation{1,Symmetric,ComplexF64,typeof(lattice_domains),typeof(nondispersive_domains),typeof(dispersive_domains),typeof(Σ)}(
		lattice_domains,
		nondispersive_domains,
		dispersive_domains,
		x,
		lattice_domain_indices,
		nondispersive_domain_indices,
		dispersive_domain_indices,
		ε,
		F,
		χ,
		laplacian,
		curlcurl,
		Σ,
		α,
		σ,
		x_half,
		σ_half,
		ω₀,
		k₁₀,
		k₂₀,
		k₃₀,
		smoothed)
end

################################################################################
# SIMULATION{1} building utilities

# used in building Simulation{1}
function _generate_half_x(lattice::Lattice{1}, x::Vector{Point{1}}, indices)
	half_x = (Vector{Point{1}}(undef,length(x)+1),)
	for i ∈ eachindex(x) half_x[1][i] = lattice[indices[i][1]-.5] end
	half_x[1][end] = lattice[indices[length(x)][1]+.5]
	return half_x
end

# 1-D hook for boundary_layer found in Simulations.jl
@inline function boundary_layer(domain::LatticeDomain{1}, x::Point{1}, dim::Integer)
	typeof(domain.boundary.shape)<:Interval || throw("1D only defined for shape=Interval, but shape=", domain.boundary.shape)
	bl1 = domain.boundary.bls[1]
	bl2 = domain.boundary.bls[2]
	σx = bl1(x)+bl2(x)
	return σx
end


################################################################################
# 1-d smoothing functions

function smooth_dielectric!(sim::Simulation{1}, num_sub_pixel::Integer = NUM_SUBPIXELS)
    indices = sim.nondispersive_domain_indices
    lattice = sim.lattice
    X = Vector{Point{1}}(undef,num_sub_pixel)
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

function smooth_pump!(sim::Simulation{1}, num_sub_pixel::Integer = NUM_SUBPIXELS)
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

function Base.getproperty(sim::Simulation{1}, sym::Symbol)
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

function Base.propertynames(sim::Simulation{1}, private=false)
	if private
		return fieldnames(Simulation)
	else
		return (:shape, :boundary, :lattice, :lattice_domain, :lattice_domains, :nondispersive_domains, :dispersive_domains, :ε, :F, :x, :χ, :ω₀, :k1₀, :k2₀, :k3₀, :lattice, propertynames(sim.lattice)...)
	end
end

################################################################################
# Pretty Printing

import ..PRINTED_COLOR_GOOD
import ..PRINTED_COLOR_WARN
import ..PRINTED_COLOR_DARK

function Base.show(io::IO,sim::Simulation{1,CLASS}) where CLASS
	if sim.smoothed[]
		printstyled(io,"Smoothed ",color=PRINTED_COLOR_GOOD)
		print(io,"1D ", CLASS)
		printstyled(io," Simulation\n",color=PRINTED_COLOR_DARK)
	else
		printstyled(io,"Unsmoothed ",color=PRINTED_COLOR_WARN)
		print(io, 1, "D ", CLASS)
		printstyled(io," Simulation\n",color=PRINTED_COLOR_DARK)
	end
    print(io,"\tn sites: ")
	printstyled(io,length(sim),"\n",color=:light_cyan)
	println(io,"\t====================")
	println(IOContext(io,:tabbed2=>true),sim.boundary)
	println(io)
	println(IOContext(io,:tabbed2=>true),sim.lattice)
	println(io)
	domains = (sim.nondispersive_domains...,sim.dispersive_domains...)
	for d ∈ eachindex(domains)
		print(io,"\t\tDomain ")
		print(io,d)
		if typeof(domains[d])<:NondispersiveDomain
			printstyled(io," nondispersive ",color=PRINTED_COLOR_DARK)
		elseif typeof(domains[d])<:DispersiveDomain
			printstyled(io," dispersive    ",color=PRINTED_COLOR_DARK)
		end
		print(io, " (",domains[d].type,"): ")
		print(io,domains[d].name)
		if typeof(domains[d])<:DispersiveDomain
			print(IOContext(io,:tabbed2=>true),"\n\t\t\t",domains[d].χ)
		end
		d < length(domains) ? println(io) : nothing
	end
end

# TODO: fix and test smoothing

function Simulation(
			ω₀::Real,
			lattice_domains::Tuple{LatticeDomain{2,Symmetric}},
			nondispersive_domains::NTuple{L,NondispersiveDomain{2}},
			dispersive_domains::NTuple{M,DispersiveDomain{2}};
			k₂₀::Real=0,
			k₃₀::Real=0,
			k₁₀::Real=sqrt(ω₀^2-k₂₀^2-k₃₀^2)
			) where {L,M}

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

	x_half = _generate_half_x(lattice_domain.lattice, x, indices)

	σ = (similar(x,ComplexF64), similar(x,ComplexF64))
	σ_half = (similar(x_half[1],ComplexF64), similar(x_half[2],ComplexF64))
	boundary_layer!(σ, lattice_domains, x, lattice_domain_indices)
	# boundary_layer!(σ_half, lattice_domains, x_half[1], vcat(lattice_domain_indices, lattice_domain_indices[end]))
	α = (1 .+ 1im*σ[1]/k₁₀, 1 .+ 1im*σ[2]/k₂₀)
	# α_half = (1 .+ 1im*σ_half[1]/k₁₀, 1 .+ 1im*σ_half[2]/k₂₀)
	return x, σ
	laplacian = Laplacian(lattice_domain.lattice,α[1],α_half[1])
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
# SIMULATION{2} building utilities

# used in building Simulation{2}
function _generate_half_x(lattice::Lattice{2}, x::Vector{Point{2,C}}, indices) where C
	half_x = (Vector{Point{2,Cartesian}}(undef,2length(x)),Vector{Point{2,Cartesian}}(undef,2length(x)))
	for i ∈ eachindex(x)
		half_x[1][i]     	   = lattice[indices[i][1]-.5, indices[i][2]   ]
		half_x[1][i+length(x)] = lattice[indices[i][1]+.5, indices[i][2]   ]
		half_x[2][i]		   = lattice[indices[i][1]   , indices[i][2]-.5]
		half_x[2][i+length(x)] = lattice[indices[i][1]   , indices[i][2]+.5]
	end
	return half_x
end

# 2-D hook for boundary_layer found in Simulations.jl
@inline function boundary_layer(domain::LatticeDomain{2}, x::Point{2}, dim::Integer)
	if typeof(domain.boundary.shape) <: AbstractQuadrilateral
		if dim==1
			bl1 = domain.boundary.bls[1]
			bl2 = domain.boundary.bls[2]
		else
			bl1 = domain.boundary.bls[3]
			bl2 = domain.boundary.bls[4]
		end
		return bl1(x)+bl2(x)
	elseif typeof(domain.boundary.shape) <: AbstractDisk
		if dim==1
			bl1 = domain.boundary.bls[1]
			bl2 = domain.boundary.bls[2]
		else
			bl1 = noBL()
			bl2 = noBL()
		end
		return bl1(x)+bl2(x)
	end
end


################################################################################
# 2-d smoothing functions

function smooth_dielectric!(sim::Simulation{2}, num_sub_pixel::Integer = NUM_SUBPIXELS)
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

function Base.propertynames(sim::Simulation{2}, private=false)
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

function Base.show(io::IO,sim::Simulation{2,CLASS}) where CLASS
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

function Simulation(
	ω₀::Real,
	domains::Vararg{Domain{1},N};
	k₂₀::Real=0,
	k₃₀::Real=0,
	k₁₀::Real=sqrt(ω₀^2-k₂₀^2-k₃₀^2)
	) where N

	lattices = map(z->z.lattice,domains)
	all(isequal.(Ref(lattices[1]),lattices)) || throw("in 1D, lattice must be shared by all domains")

	starts = map(d->d.boundary.shape.start,domains)
	stops =  map(d->d.boundary.shape.stop ,domains)
	start = min(starts...,stops...)
	stop  = max(starts...,stops...)
	a::Float64 = stop-start
	nx = ceil(Int,a/domains[1].lattice.dx)
	dx = a/nx

	lat = Lattice(dx;origin=start+dx/2)

	domains = map(d->Domain(d.boundary,lat,d.dielectric,d.pump,d.χ;type=d.type,name=d.name),domains)

	smoothed = Ref(false)

	domain_index = Int[]
	indices = CartesianIndex{1}[]
	nnm = (Int[],)
	nnp = (Int[],)
	name = Symbol[]
	x = Point{1}[]
	ε = ComplexF64[]
	F = Float64[]
	interior = BitArray([])
	exterior = falses(length(x))
	surface = BitArray([])
	bulk = BitArray([])
	corner = BitArray([])

	append_domains!(domain_index,nnm,nnp,name,x,indices,ε,F,interior,bulk,surface,corner,domains)

	removed = get_sites_to_remove(domains,domain_index,x)

	trim_nn!(nnm,nnp,removed,bulk)

	remove_sites!((domain_index,name,x,indices,ε,F,interior,bulk,surface,corner,nnm,nnp),removed)

	x_half = generate_half_x(domains,domain_index,x,indices)

	domain_wall = generate_and_domain_wall!(interior,bulk,surface,corner,domains,domain_index,x,indices)

	σ = boundary_layer(domains,x,domain_index)
	σ_half = boundary_layer(domains,x_half[1],vcat(domain_index,domain_index[end]))
	α = 1 .+ 1im*σ[1]/k₁₀
	α_half = 1 .+ 1im*σ_half[1]/k₁₀

	curlcurl = Curlcurl(domains[1].lattice,α,α_half,nnm,nnp,indices,interior,surface,domain_wall)
	Σ = SelfEnergy(domains,surface,domain_index,α_half,a,nnm,nnp,indices,interior,domain_wall)

	return Simulation{1,Symmetric,ComplexF64,typeof(domains),2,typeof(Σ.f)}(
		smoothed,
		ε,
		F,
		curlcurl,
		Σ,
		domains,
		domain_index,
		x,
		x_half,
		(α,),
		σ,
		σ_half,
		name,
		ω₀,
		k₁₀,
		k₂₀,
		k₃₀)
end


function generate_and_domain_wall!(interior,bulk,surface,corner,domains,domain_index,x::Array{Point{1}},indices)
	lattices = map(z->z.lattice,domains)
	shapes = map(z->z.boundary.shape,domains)
	DX = map(z->z.dx,lattices)
	domain_wall = falses(length(x))
	lattice_wall = falses(length(x))
	for i ∈ eachindex(x)
		if interior[i]
			ds = domain_index[i]
			xt = lattices[ds][indices[i]-CartesianIndex(1,)]
			dom1 = which_domain(domains,shapes,xt)
			xt = lattices[ds][indices[i]+CartesianIndex(1,)]
			dom2 = which_domain(domains,shapes,xt)

			domain_wall[i] = !(domain_index[i]==dom1==dom2)
			domain_wall[i] = domain_wall[i] && (dom1!==0 && dom2!==0)
			if domain_wall[i]
				corner[i] = false
				surface[i] = false
				bulk[i] = true
			end
		end
	end
	return domain_wall
end

function generate_half_x(domains,domain_index,x::Array{Point{1}},indices)
	half_x = (Array{Point{1},1}(undef,length(x)+1),)
	lattices = map(z->z.lattice,domains)
	for i ∈ 1:length(x)
		ds = domain_index[i]
		half_x[1][i] = lattices[ds][indices[i][1]-.5]
		if i==length(x)
			half_x[1][i+1] = lattices[ds][indices[i][1]+.5]
		end
	end
	return half_x
end

function smooth!(sim::Simulation{1},num_sub_pixel::Int=NUM_SUB_PIXEL)
    domains = sim.domains
    domain_wall = sim.domain_wall
	lattices = map(z->z.lattice,sim.domains)
    X = Array{Point{1}}(undef,num_sub_pixel)
    E = Array{ComplexF64}(undef,num_sub_pixel)
    F = Array{Float64}(undef,num_sub_pixel)
	dw_inds = findall(domain_wall)
    for dw ∈ eachindex(dw_inds)
		d = dw_inds[dw]
		ds = sim.domain_index[d]
	    I = Tuple(sim.indices[d])
        is = LinRange(I-.5,I+.5,num_sub_pixel)
        for i ∈ eachindex(is), j ∈ eachindex(js)
            X[i] = lattices[ds][is[i]]
        end
        simulation_dielectric!(E,domains,X)
        simulation_pump!(F,domains,X)
        sim.ε[d] = mean(E)
        sim.F[d] = mean(F)
    end
	sim.smoothed[] = true
    return nothing
end

function boundary_layer(domain,x::Point{1},dim::Int)
	typeof(domain.boundary.shape)<:Interval || throw("1D only defined for shape=Interval, but shape=",domain.boundary.shape)
	N = 2
	bl1 = domain.boundary.bls[1]
	bl2 = domain.boundary.bls[2]
	σx = bl1(x)+bl2(x)
	return σx
end

function Base.conj(sim::Simulation)
	return Simulation(sim.ω₀,map(conj,sim.domains)...)
end

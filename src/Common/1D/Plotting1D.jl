import ..Common.SHAPE_FILL_ALPHA
import ..Common.SHAPE_COLOR

# Simulation
@recipe function f(sim::Simulation{1}, subplots=[1,2,3])
    x = map(a->a.vec[1],sim.x)
    n = sqrt.(sim.ε)

    layout --> (3,1)
    legend --> false

	for i ∈ 1:2
		if !(typeof(sim.boundary.bls[i])<:noBL)
			for j ∈  subplots
				@series begin
					subplot := j
					seriestype --> :path
					alpha --> 0
					if typeof(sim.boundary.bls[i])<:PML
						fillcolor --> :blue
					elseif typeof(sim.boundary.bls[i])<:cPML
						fillcolor --> :red
					end
					fillrange --> 0
					fillalpha --> SHAPE_FILL_ALPHA
					([sim.boundary.bls[i].start,
					  sim.boundary.bls[i].stop,
					  sim.boundary.bls[i].stop,
					  sim.boundary.bls[i].start,
					  sim.boundary.bls[i].start],

					  [-1000, -1000, 1000, 1000, -1000])
				end
			end
		end
	end
    @series begin
		title --> L"{\rm Scattering\ Structure}"
		ylabel --> L"{\rm real}(n)"
        subplot := subplots[1]
		ylims --> (min(1,minimum(real(n))), max(1,maximum(real(n)))) .+ (-.1,.1)
        x, real(n)
    end
    @series begin
		ylabel --> L"{\rm imag}(n)"
        subplot := subplots[2]
		ylims --> (-1,1).*max(1,maximum(abs.(imag(n)))) .+ (-.1,.1)
        x, imag(n)
    end
    @series begin
		xlabel --> L"x"
		ylabel --> L"{\rm pump\ F}"
		subplot := subplots[3]
		ylims --> (-1,1).*max(1,maximum(abs.(sim.F))) .+ (-.1,.1)
        x, sim.F
    end
end


# Scalar Field
@recipe function f(e::ScalarField{1}, by::Function)
    x = map(a->a.vec[1],e.pos)

	upper_lim = maximum(abs∘by,e.val)
	if by ∈ (abs,abs2)
		lower_lim = zero(upper_lim)
	else
		lower_lim = -upper_lim
	end

	m = size(e,2)
    layout --> (m,1)
    legend --> false
	ylims --> (lower_lim,upper_lim)
	for μ ∈ 1:m
	    @series begin
			subplot --> μ
			m>1 ? title --> latexstring("{\\rm Field}\\ $μ") : nothing
	        x, by.(e(μ).val)
	    end
	end
end


@recipe function f(sim::Simulation{1}, e::ScalarField{1})
	legend --> false
	m = size(e,2)
	layout --> (3,1+m)

	@series begin
		sim, [1,1+1+m,1+2(1+m)]
	end
    for μ ∈ 1:m
		@series begin
			subplot --> μ+1
			e(μ), real
		end
		@series begin
			subplot --> μ+1 + m+1
			e(μ), imag
		end
		@series begin
			subplot --> μ+1 + 2(m+1)
			e(μ), abs2
		end
	end
end

# Electric Field
@recipe function f(e::ElectricField{1}, by::Function)
    x = map(a->a.vec[1],e.pos)
    perm = sortperm(x)
	x = x[perm]

	upper_lim = maximum(abs∘by,e.val)
	if by ∈ (abs,abs2)
		lower_lim = zero(upper_lim)
	else
		lower_lim = -upper_lim
	end

	m = size(e,2)
    layout --> (3,m)
    legend --> false
	ylims --> (lower_lim,upper_lim)
	for μ ∈ 1:m
	    @series begin
			title --> latexstring("{\\rm Field}\\ $μ")
	        ylabel --> L"E_x"
	        subplot := 0m+μ
	        x, by.(e(μ).x[perm])
	    end
	    @series begin
	        ylabel --> L"E_y"
	        subplot := 1m+μ
	        x, by.(e(μ).y[perm])
	    end
	    @series begin
			xlabel --> L"x"
	        ylabel --> L"E_z"
	        subplot := 2m+μ
	        x, by.(e(μ).z[perm])
	    end
	end
end


@recipe function f(sim::Simulation{1},e::ElectricField{1}, by::Function)
    x = map(a->a.vec[1],e.pos)
    perm = sortperm(x)
    x = x[perm]

	legend --> false

	m = 1+size(e,2)
	layout --> (3,m)

	n = sqrt.(sim.ε[perm])
	@series begin
		title --> L"{\rm Scattering\ Structure}"
		ylabel --> L"{\rm real}(n)"
        subplot := 0m+1
        x, real(n)
    end
    @series begin
		ylabel --> L"{\rm imag}(n)"
        subplot := 1m+1
        x, imag(n)
    end
    @series begin
		xlabel --> L"x"
		ylabel --> L"{\rm pump\ F}"
		subplot := 2m+1
        x, sim.F[perm]
    end

	upper_lim = maximum(abs∘by,e.val)
	if by ∈ (abs,abs2)
		lower_lim = zero(upper_lim)
	else
		lower_lim = -upper_lim
	end
	ylims --> (lower_lim,upper_lim)

	for μ ∈ 1:m-1
	    @series begin
			title --> latexstring("{\\rm Field}\\ $μ")
	        ylabel --> L"E_x"
	        subplot := 0m+μ+1
	        x, by.(e(μ).x[perm])
	    end
	    @series begin
	        ylabel --> L"E_y"
	        subplot := 1m+μ+1
	        x, by.(e(μ).y[perm])
	    end
	    @series begin
			xlabel --> L"x"
	        ylabel --> L"E_z"
	        subplot := 2m+μ+1
	        x, by.(e(μ).z[perm])
	    end

	end
end

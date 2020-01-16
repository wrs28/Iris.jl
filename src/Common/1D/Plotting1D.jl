# Simulation
@recipe function f(sim::Simulation{1})
    x = map(a->a.vec[1],sim.x)
    perm = sortperm(x)
    x = x[perm]
    n = sqrt.(sim.ε[perm])

    layout --> (3,1)
    legend --> false

    @series begin
		title --> L"{\rm Scattering\ Structure}"
		ylabel --> L"{\rm real}(n)"
        subplot := 1
        x, real(n)
    end
    @series begin
		ylabel --> L"{\rm imag}(n)"
        subplot := 2
        x, imag(n)
    end
    @series begin
		xlabel --> L"x"
		ylabel --> L"{\rm pump\ F}"
		subplot := 3
        x, sim.F[perm]
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



# Scalar Field
@recipe function f(e::ScalarField{1}, by::Function)
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
    layout --> (m,1)
    legend --> false
	ylims --> (lower_lim,upper_lim)
	for μ ∈ 1:m
	    @series begin
			title --> latexstring("{\\rm Field}\\ $μ")
	        ylabel --> L"\Psi"
	        subplot := μ
	        x, by.(e(μ).val[perm])
	    end
	end
end


@recipe function f(sim::Simulation{1},e::ScalarField{1}, by::Function)
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
	end
end

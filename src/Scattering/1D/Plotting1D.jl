using LaTeXStrings
using RecipesBase

# Helmholtz scattering solution only
@recipe function f(sct::ScatteringSolution{1,1}, by::Function)

	upper_lim = maximum(abs∘by,sct.tot.val)
	upper_lim = max(upper_lim,maximum(abs∘by,sct.sct.val))
	upper_lim = max(upper_lim,maximum(abs∘by,sct.inc.val))
	if by ∈ (abs,abs2)
		lower_lim = zero(upper_lim)
	else
		lower_lim = -upper_lim
	end

	layout --> (3,1)
	legend --> false
	ylims --> (lower_lim,upper_lim)

	@series begin
		title := L"{\rm Incident\ Field}"
		ylabel := L"E_x"
		subplot := 1
		sct.incident, by
	end
	@series begin
		title := L"{\rm Scattered\ Field}"
		subplot := 2
		sct.scattered, by
	end
	@series begin
		title := L"{\rm Total\ Field}"
		subplot := 3
		sct.total, by
	end
end

# Helmotz scattering + simulation
@recipe function f(sim::Simulation{1},sct::ScatteringSolution{1,1}, by::Function)

	upper_lim = maximum(abs∘by,sct.tot.val)
	upper_lim = max(upper_lim,maximum(abs∘by,sct.sct.val))
	upper_lim = max(upper_lim,maximum(abs∘by,sct.inc.val))
	if by ∈ (abs,abs2)
		lower_lim = zero(upper_lim)
	else
		lower_lim = -upper_lim
	end

	layout --> (3,2)
	legend --> false

	@series begin
		sim, [1,3,5]
	end
	ylims --> (lower_lim,upper_lim)
	@series begin
		title := L"{\rm Incident\ Field}"
		ylabel := L"E_x"
		subplot := 2
		sct.incident, by
	end
	@series begin
		title := L"{\rm Scattered\ Field}"
		subplot := 4
		sct.scattered, by
	end
	@series begin
		title := L"{\rm Total\ Field}"
		subplot := 6
		sct.total, by
	end
end





@recipe function f(sct::ScatteringSolution{1,3},by::Function)
    x = map(a->a.vec[1],sct.tot.pos)
    perm = sortperm(x)
    x = x[perm]

	upper_lim = maximum(abs∘by,sct.tot.val)
	upper_lim = max(upper_lim,maximum(abs∘by,sct.sct.val))
	upper_lim = max(upper_lim,maximum(abs∘by,sct.inc.val))
    if by ∈ (abs,abs2)
        lower_lim = zero(upper_lim)
    else
        lower_lim = -upper_lim
    end

    layout --> (3,3)
    legend --> false
    ylims --> (lower_lim,upper_lim)

    @series begin
        title := L"{\rm Incident\ Field}"
        ylabel := L"E_x"
        subplot := 1
        x, by.(sct.inc.x[perm])
    end
    @series begin
        ylabel := L"E_y"
        subplot := 4
        x, by.(sct.inc.y[perm])
    end
    @series begin
        xlabel := L"x"
        ylabel := L"E_z"
        subplot := 7
        x, by.(sct.inc.z[perm])
    end

    @series begin
        title := L"{\rm Scattered\ Field}"
        subplot := 2
        x, by.(sct.sct.x[perm])
    end
    @series begin
        subplot := 5
        x, by.(sct.sct.y[perm])
    end
    @series begin
        xlabel := L"x"
        subplot := 8
        x, by.(sct.sct.z[perm])
    end

    @series begin
        title := L"{\rm Total\ Field}"
        subplot := 3
        x, by.(sct.tot.x[perm])
    end
    @series begin
        subplot := 6
        x, by.(sct.tot.y[perm])
    end
    @series begin
        xlabel := L"x"
        subplot := 9
        x, by.(sct.tot.z[perm])
    end
end

# Scattering + Simulation
@recipe function f(sim::Simulation{1},sct::ScatteringSolution{1},by::Function)
    x = map(a->a.vec[1],sct.tot.pos)
    perm = sortperm(x)
    x = x[perm]

    layout --> (3,4)
    legend --> false

    n = sqrt.(sim.ε[perm])
    @series begin
		title --> L"{\rm Scattering\ Structure}"
		ylabel --> L"{\rm real}(n)"
        subplot := 1
        x, real(n)
    end
    @series begin
		ylabel --> L"{\rm imag}(n)"
        subplot := 5
        x, imag(n)
    end
    @series begin
		xlabel --> L"x"
		ylabel --> L"{\rm pump\ F}"
		subplot := 9
        x, sim.F[perm]
    end

    upper_lim = maximum(abs∘by,sct.tot.val)
	upper_lim = max(upper_lim,maximum(abs∘by,sct.sct.val))
	upper_lim = max(upper_lim,maximum(abs∘by,sct.inc.val))
    if by ∈ (abs,abs2)
        lower_lim = zero(upper_lim)
    else
        lower_lim = -upper_lim
    end

    ylims --> (lower_lim,upper_lim)

    @series begin
        title := L"{\rm Incident\ Field}"
        ylabel := L"E_x"
        subplot := 2
        x, by.(sct.inc.x[perm])
    end
    @series begin
        ylabel := L"E_y"
        subplot := 6
        x, by.(sct.inc.y[perm])
    end
    @series begin
        xlabel := L"x"
        ylabel := L"E_z"
        subplot := 10
        x, by.(sct.inc.z[perm])
    end

    @series begin
        title := L"{\rm Scattered\ Field}"
		ylabel := L"E_x"
        subplot := 3
        x, by.(sct.sct.x[perm])
    end
    @series begin
		ylabel := L"E_y"
        subplot := 7
        x, by.(sct.sct.y[perm])
    end
    @series begin
		xlabel := L"x"
		ylabel := L"E_z"
        subplot := 11
        x, by.(sct.sct.z[perm])
    end

    @series begin
        title := L"{\rm Total\ Field}"
		ylabel := L"E_x"
        subplot := 4
        x, by.(sct.tot.x[perm])
    end
    @series begin
		ylabel := L"E_y"
        subplot := 8
        x, by.(sct.tot.y[perm])
    end
    @series begin
        xlabel := L"x"
		ylabel := L"E_z"
        subplot := 12
        x, by.(sct.tot.z[perm])
    end
end

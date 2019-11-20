using LaTeXStrings
using RecipesBase

@recipe f(mls::Maxwell_LS{1};by=abs2) = mls.solution,by
@recipe f(mls::Maxwell_LS{1},sim::Simulation{1};by=abs2) = sim,mls.solution,by
@recipe f(sim::Simulation{1},mls::Maxwell_LS{1};by=abs2) = sim,mls.solution,by

@recipe f(nls::Maxwell_NLS{1};by=abs2) = nls.solution,by
@recipe f(nls::Maxwell_NLS{1},sim::Simulation{1};by=abs2) = sim,nls.solution,by
@recipe f(sim::Simulation{1},nls::Maxwell_NLS{1};by=abs2) = sim,nls.solution,by

# Scattering
@recipe f(sct::ScatteringSolution;by=abs2) = sct,by
@recipe f(by::Function, sct::ScatteringSolution) = sct,by

@recipe function f(sct::ScatteringSolution{1},by::Function)
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
@recipe f(sim::Simulation{1},sct::ScatteringSolution{1};by=abs2) = sim,sct,by
@recipe f(sct::ScatteringSolution{1},sim::Simulation{1};by=abs2) = sim,sct,by
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

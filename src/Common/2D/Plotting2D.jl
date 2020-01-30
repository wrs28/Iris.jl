import ..Common.MARKERSIZE_SCALE
import ..Common.MARKERSHAPE
import ..Symmetric

# Cartesian Symmetric Simulation
@recipe function f(sim::Simulation{2,Symmetric,C,Tuple{LD}}) where {C,LD<:LatticeDomain{2,Symmetric,Cartesian}}
	lattice = sim.lattice_domain.lattice
	indices = sim.lattice_domain.indices
	imin, imax = extrema(map(ld->ld[1],indices))
	jmin, jmax = extrema(map(ld->ld[2],indices))
	x = [lattice[i,0].x for i ∈ imin:imax]
	y = [lattice[0,j].y for j ∈ jmin:jmax]

	n = permutedims(reshape(sqrt.(sim.ε[:]),length(x),length(y)))

    layout --> (1,3)
    legend --> false

	climr = maximum(real(n))
	climi = max(maximum(abs.(imag(n))),1e-5)
	climf = max(maximum(sim.F),1e-5)

    @series begin
		title --> L"{\rm Real}(n)"
        subplot := 1
        seriestype --> :heatmap
		seriescolor --> :sequential
		aspectratio --> 1
		clims --> (1,climr)
        x, y, real(n)
    end
	@series begin
		title --> L"{\rm Imag}(n)"
		subplot := 2
        seriestype --> :heatmap
		seriescolor --> :diverging
		aspectratio --> 1
		clims --> (-climi,climi)
        x, y, imag(n)
    end
	@series begin
		title --> L"{\rm Pump\ }F"
        subplot := 3
        seriestype --> :heatmap
		seriescolor --> :diverging
		aspectratio --> 1
		subplot := 3
		zcolor --> sim.F[:]
		clims --> (-climf,climf)
        y, x, reshape(sim.F[:],length(x),:)
    end
end

# Generic Simulation
@recipe function f(sim::Simulation{2})
    n = sqrt.(sim.ε[:])

    layout := (1,3)
    legend --> false

    @series begin
		title --> L"{\rm Real}(n)"
        subplot := 1
        seriestype --> :scatter
        markersize --> MARKERSIZE_SCALE/sqrt(length(sim.x))
        markerstrokealpha --> 0
        shape --> MARKERSHAPE
		seriescolor --> :sequential
		aspectratio --> 1
		zcolor --> real(n)
		clims --> (1,3)
        sim.x
    end
	@series begin
		title --> L"{\rm Imag}(n)"
        subplot := 2
        seriestype --> :scatter
        markersize --> MARKERSIZE_SCALE/sqrt(length(sim.x))
        markerstrokealpha --> 0
        shape --> MARKERSHAPE
		seriescolor --> :diverging
		aspectratio --> 1
		zcolor --> imag(n)
		clims --> (-1,1)
        sim.x
    end
	@series begin
		title --> L"{\rm Pump\ }F"
        subplot := 3
        seriestype --> :scatter
        markersize --> MARKERSIZE_SCALE/sqrt(length(sim.x))
        markerstrokealpha --> 0
        shape --> MARKERSHAPE
		seriescolor --> :diverging
		aspectratio --> 1
		subplot := 3
		zcolor --> sim.F[:]
		clims --> (-1,1)
        sim.x
    end
end

# Cartesian Symmetric ScalarField
@recipe function f(sim::Simulation{2,Symmetric,C,Tuple{LD}}, ψ::ScalarField{2}) where {C,LD<:LatticeDomain{2,Symmetric,Cartesian}}
	m = size(ψ,2)
	layout --> (1+m,3)
    legend --> false
	@series begin
		sim
	end
	for μ ∈ 1:m
		ψμ = ψ(μ)
		clim = max(maximum(abs.(real(ψμ))), maximum(abs.(real(ψμ))))
	    @series begin
			title --> L"{\rm Real}(\psi)"
	        subplot := 3μ+1
			sim, ψμ, real
	    end
		@series begin
			title --> L"{\rm Imag}(\psi)"
	        subplot := 3μ+2
	        sim, ψμ, imag
	    end
		@series begin
			title --> L"|\psi|^2"
	        subplot := 3μ+3
	        sim, ψμ, abs2
	    end
	end
end

# Cartesian Symmetric ScalarField by function
@recipe function f(sim::Simulation{2,Symmetric,C,Tuple{LD}}, ψ::ScalarField{2}, by::Function) where {C,LD<:LatticeDomain{2,Symmetric,Cartesian}}

	lattice = sim.lattice_domain.lattice
	indices = sim.lattice_domain.indices
	imin, imax = extrema(map(ld->ld[1],indices))
	jmin, jmax = extrema(map(ld->ld[2],indices))
	x = [lattice[i,0].x for i ∈ imin:imax]
	y = [lattice[0,j].y for j ∈ jmin:jmax]

	legend --> false

	for μ ∈ 1:size(ψ,2)
		ψμ = ψ(μ)
		clim = max(maximum(abs.(by.(ψμ))), maximum(abs.(by.(ψμ))))
	    @series begin
	        seriestype --> :heatmap
			aspectratio --> 1
			if by ∈ (abs, abs2)
				clims --> (0,clim)
				seriescolor --> :sequential
			else
				clims --> (-clim,clim)
				seriescolor --> :diverging
			end
	        x, y, by.(ψμ)
	    end
	end
end

# Generic ScalarField
@recipe function f(ψ::ScalarField{2})
	m = size(ψ,2)

	layout := (m,3)
    legend --> false
	for μ ∈ 1:m
		ψμ = ψ(μ)

		clim = max(maximum(abs.(real(ψμ))), maximum(abs.(real(ψμ))))
	    @series begin
			title --> L"{\rm Real}(\psi)"
	        subplot := 3*(μ-1)+1
	        seriestype --> :scatter
	        markersize --> MARKERSIZE_SCALE/sqrt(length(ψ.positions))
	        markerstrokealpha --> 0
	        shape --> MARKERSHAPE
			seriescolor --> :diverging
			aspectratio --> 1
			zcolor --> real(ψμ)
			clims --> (-clim,clim)
	        ψ.positions
	    end
		@series begin
			title --> L"{\rm Imag}(\psi)"
	        subplot := 3*(μ-1)+2
	        seriestype --> :scatter
	        markersize --> MARKERSIZE_SCALE/sqrt(length(ψ.positions))
	        markerstrokealpha --> 0
	        shape --> MARKERSHAPE
			seriescolor --> :diverging
			aspectratio --> 1
			zcolor --> imag(ψμ)
			clims --> (-clim,clim)
	        ψ.positions
	    end
		@series begin
			title --> L"|\psi|^2"
	        subplot := 3*(μ-1)+3
	        seriestype --> :scatter
	        markersize --> MARKERSIZE_SCALE/sqrt(length(ψ.positions))
	        markerstrokealpha --> 0
	        shape --> MARKERSHAPE
			seriescolor --> :sequential
			aspectratio --> 1
			zcolor --> abs2.(ψμ)
			clims --> (0, maximum(abs2.(ψμ)))
	        ψ.positions
	    end
	end
end


# Electric Field
# @recipe function f(e::ElectricField{1}, by::Function)
    # x = map(a->a.vec[1],e.pos)
    # perm = sortperm(x)
	# x = x[perm]

	# upper_lim = maximum(abs∘by,e.val)
	# if by ∈ (abs,abs2)
		# lower_lim = zero(upper_lim)
	# else
		# lower_lim = -upper_lim
	# end
#
# 	m = size(e,2)
#     layout --> (3,m)
#     legend --> false
# 	ylims --> (lower_lim,upper_lim)
# 	for μ ∈ 1:m
# 	    @series begin
# 			title --> latexstring("{\\rm Field}\\ $μ")
# 	        ylabel --> L"E_x"
# 	        subplot := 0m+μ
# 	        x, by.(e(μ).x[perm])
# 	    end
# 	    @series begin
# 	        ylabel --> L"E_y"
# 	        subplot := 1m+μ
# 	        x, by.(e(μ).y[perm])
# 	    end
# 	    @series begin
# 			xlabel --> L"x"
# 	        ylabel --> L"E_z"
# 	        subplot := 2m+μ
# 	        x, by.(e(μ).z[perm])
# 	    end
# 	end
# end
#
#
# @recipe function f(sim::Simulation{1},e::ElectricField{1}, by::Function)
#     x = map(a->a.vec[1],e.pos)
#     perm = sortperm(x)
#     x = x[perm]
#
# 	legend --> false
#
# 	m = 1+size(e,2)
# 	layout --> (3,m)
#
# 	n = sqrt.(sim.ε[perm])
# 	@series begin
# 		title --> L"{\rm Scattering\ Structure}"
# 		ylabel --> L"{\rm real}(n)"
#         subplot := 0m+1
#         x, real(n)
#     end
#     @series begin
# 		ylabel --> L"{\rm imag}(n)"
#         subplot := 1m+1
#         x, imag(n)
#     end
#     @series begin
# 		xlabel --> L"x"
# 		ylabel --> L"{\rm pump\ F}"
# 		subplot := 2m+1
#         x, sim.F[perm]
#     end
#
# 	upper_lim = maximum(abs∘by,e.val)
# 	if by ∈ (abs,abs2)
# 		lower_lim = zero(upper_lim)
# 	else
# 		lower_lim = -upper_lim
# 	end
# 	ylims --> (lower_lim,upper_lim)
#
# 	for μ ∈ 1:m-1
# 	    @series begin
# 			title --> latexstring("{\\rm Field}\\ $μ")
# 	        ylabel --> L"E_x"
# 	        subplot := 0m+μ+1
# 	        x, by.(e(μ).x[perm])
# 	    end
# 	    @series begin
# 	        ylabel --> L"E_y"
# 	        subplot := 1m+μ+1
# 	        x, by.(e(μ).y[perm])
# 	    end
# 	    @series begin
# 			xlabel --> L"x"
# 	        ylabel --> L"E_z"
# 	        subplot := 2m+μ+1
# 	        x, by.(e(μ).z[perm])
# 	    end
#
# 	end
# end
